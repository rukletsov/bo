
/******************************************************************************

    d25_active_contours.hpp, v 1.0.0 2011.03.17

	Modification of the approach by Ye Duan and Hong Qin, "2.5D Active Contour
	for Surface Reconstruction", Proceedings of the 8th Fall Workshop on Vision,
	Modeling and Visualization (VMV 2003), Munich, Germany, November 19-21, 2003,
	pages 431 -- 439.

    Copyright (c) 2009-2011, Dzmitry Hlindzich <hlindzich@gmail.com>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
	    notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
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

#include "point.hpp"
#include "triangle.hpp"
#include "mesh.hpp"

namespace methods {

namespace {

/*! \class HVertexContainer
	\brief A general container of vertices
*/
class HVertexContainer;

} //anonymous namespace

namespace surfaces {

/*! \struct HPointSeed
	\brief An elementary point item
*/
struct HPointSeed
{	
	/*! 3D Vertex */
	common::Point3<float> p;

	/*! Visits flag */
	bool isVisited;

	/*! Nodes flag */
	bool isNode;

	/*! Access operator */
	inline float operator [] (const int t) const 
	{
		return t==0?p.x:(t==1?p.y:p.z);
	}

	/*! Comparison operator */
	inline bool operator == (const HPointSeed &other) const
	{
		return (p==other.p)&&(isNode==other.isNode)&&(isVisited==other.isVisited);
	}
};


/*! \class HEdgeSeed
	\brief An elementary edge propagating item
*/
struct  HEdgeSeed
{
	/*! Default constructor */
	HEdgeSeed():p1(0),p2(0){}

	/*! Constructor */
	HEdgeSeed(HPointSeed* p1, HPointSeed* p2):p1(p1), p2(p2){}

	/*! First node of the edge */
	HPointSeed* p1;

	/*! Second node of the edge */
	HPointSeed* p2;

	/*! Vector of the propagation direction */
	common::Point3<float> propagationVector;

	/*! Comparison operator */
	bool operator == (const HEdgeSeed& other) const
	{
		return  ((p1==other.p1&&p2==other.p2)||(p1==other.p2&&p2==other.p1))&&
			(propagationVector==other.propagationVector);
	}

	/*! Swaps the vertices of the edge */ 
	void swap()
	{
		HPointSeed* tmp=p1;
		p1=p2;
		p2=tmp;
	}
};


/*! \class HTriangleSeed
	\brief An elementary surface item (mesh polygon)
*/
struct  HTriangleSeed
{
	/*! First triangle's vertex */
	HPointSeed* p1;

	/*! Second triangle's vertex */
	HPointSeed* p2;

	/*! Third triangle's vertex */
	HPointSeed* p3;
};


/*! \class D25ActiveContours
	\brief 2.5D active contour based mesh reconstruction
	\author Dzmitry Hlindzich
	\date 2009-2011
	\details Modification of the approach proposed by Ye Duan and Hong Qin,
	"2.5D Active Contour for Surface Reconstruction", Proceedings of the 8th
	Fall Workshop on Vision, Modeling and Visualization (VMV 2003), Munich,
	Germany, November 19-21, 2003, pages 431 -- 439.
*/
class  D25ActiveContours
{
public:

	/*! A constructor */
	D25ActiveContours();

	/*! A destructor */
	~D25ActiveContours();

	/*! Initialization of all the parameters based on the minimal length of the initial triangle's side
		\param minInitDistance The minimal allowed length of the initial triangle's side
		\see minInitDistance
	*/
	void initMeshScale(float minInitDistance);

	/*! Sets the minimal allowed length of the initial triangle's side
		\param minInitDistance The minimal allowed length of the initial triangle's side
		\see minInitDistance
	*/
	void setMinInitDistance(float minInitDistance);

	/*! Sets the maximal allowed length of the initial triangle's side
		\param maxInitDistance The maximal allowed length of the initial triangle's side
		\see maxInitDistance
	*/
	void setMaxInitDistance(float maxInitDistance);

	/*! Sets maximum allowed distance between a propagated point to the correspondent projection vertex in 3D
		\param maxProjectionNodeDistance The maximal allowed distance between a point and it's node-projection in 3D
		\see maxProjectionNodeDistance
	*/
	void setMaxProjectionNodeDistance(float maxProjectionNodeDistance);

	/*! Sets the radius of the neighborhood system used for PCA-based normal vector calculation
		\param normalNeighborhoodRadius The radius of the neighborhood system used for PCA-baesd normal vector calculation
		\see normalNeighborhoodRadius
	*/
	void setNormalNeighborhoodRadius(float normalNeighborhoodRadius);

	/*! Sets the maximum analyzed depth of the points cloud layer. Thicker layers will be interpreted as several surfaces 
		\param maxSurfaceDepth The maximum analyzed depth of the points cloud layer
		\see maxSurfaceDepth
	*/
	void setMaxSurfaceDepth(float maxSurfaceDepth);

	/*! Sets the maximum value of angles between mesh edges that will be excluded during the mesh construction. Reduces the number of "thin" polygons in the mesh. Takes the values from 0.0 to 1.0 (angles in cos values, e.g. 0 degrees -> 1.0, 90 degress -> 0.0 )
		\param maxExcludedAngle The maximum value of angles between mesh edges that will be excluded during the mesh construction
		\see maxExcludedAngle
	*/
	void setMaxExcludedAngle(float maxExcludedAngle);

	/*! Sets the linear proportion between PCA-normal and triangle-based normal in the surface normal calculation procedure. Values from 0.0 to 1.0 only. 0.0 for the triangle-based normals, 1.0 for the PCA-based ones.
		\param faceSurfaceFactor The linear proportion between PCA-normal and triangle-based normals
	*/
	void setFaceSurfaceFactor(float faceSurfaceFactor);

	/*! Sets the maximum cos-value of the angles between adjacent frozen edges, that allows post-stitching of the edges
		\param maxStitchedAngle the maximum cos-value of the angles between adjacent frozen edges, that allows post-stitching of the edges
		\see maxStitchedAngle
	*/
	void setMaxStitchedAngle(float maxStitchedAngle);

	/*! Vertices vector getter
		\return Vertices vector
	*/
	const std::vector<HPointSeed>* getVerticesVector();

	/*! Active edges list getter
		\return Active edges list
		\see activeEdges
	*/
	const std::list<HEdgeSeed>* getActiveEdgeList();
	/*! Frozen edges list getter
		\return Frozen edges list
		\see frozenEdges
	*/
	const std::list<HEdgeSeed>* getFrozenEdgeList();
	/*! Triangles list getter
		\return Triangles list
		\see activeEdges
	*/
	const std::list<HTriangleSeed>* getTriangleList();

	/*! The main mesh construction procedure.
		\param vertexList The list of vertices in 3D
		\return Reconstructed mesh
	*/
	common::Mesh buildMesh(std::list<common::Point3<float>> &vertexList);

	/*! Loads points cloud into the internal vertices structure
		\param vertexList The list of vertices in 3D
	*/
	void loadVertices(std::list<common::Point3<float>> &vertexList);
	
	/*! Performs one growing iteration. Returns true if further growing is possible otherwise returns false
		\return true if further growing is possible otherwise returns false
	*/
	bool growStep();


protected:

	/*!Generalizes an init triangle seed from unvisited elements of \p vertices if possible.  Visits the points
	*/
	void modelInit();
	
	/*! Tries to propagate an element of the active contour. The surface can be propagated if \p activeEdges is not empty
	*/
	void modelGrow();
	
	/*! Performs "on-the-fly" stitching (adding one more triangle) of the given edge \p e to the mesh \p triangles if the angle between them is small enough
		\param e The stitched edge
	*/
	void edgeStitch(HEdgeSeed e);

	/*! Performs one step of "post-stitching" (connection of adjacent frozen edges into triangles if the angle between them is less than \p maxStitchedAngle) with possible addition of new active edges and deletion of the stitched frozen ones
	*/
	void postStitch();


	/*! Finds the closest point to the given \p ps that lies in the distance interval from \p minInitDistance to \p maxInitDistance. Searches among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true 
		\param ps The origin vertex
		\param checkNodes A flag. Search among the nodes if true
		\param checkVisited A flag. Search among the visited vertices if true
		\return A pointer to the resulting vertex
		\see maxInitDistance
		\see minInitDistance
	*/
	HPointSeed* getClosestPoint(const HPointSeed &ps, bool checkNodes, bool checkVisited);
	
	/*! Finds the closest point P to the given \p ps such that P, \p ps1 and \p ps2 are non-collinear. Searches among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true 
		\param ps The origin vertex
		\param ps1 First reference vertex
		\param ps2 Second reference vertex
		\param checkNodes A flag. Search among the nodes if true
		\param checkVisited A flag. Search among the visited vertices if true
		\return A pointer to the resulting vertex
	*/
	HPointSeed* getClosestNoncollinearPoint(const HPointSeed &ps, const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited);
	
	/*! Finds point P that minimizes F(P, \p ps1, \p ps2)=abs(|ps1-ps2|-|P-ps2|)+abs(|ps1-ps2|-|P-ps1|), |P-ps1|,|P-ps2|<\p maxInitDistance. Searches among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true 
		\param ps1 First reference vertex
		\param ps2 Second reference vertex
		\param checkNodes A flag. Search among the nodes if true
		\param checkVisited A flag. Search among the visited vertices if true
		\return A pointer to the resulting vertex
		\see maxInitDistance
	*/
	HPointSeed* getClosestMinFuncPoint(const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited);
	
	/*! Calculates the Euclidean distance between \p ps1 and \p ps2
		\param ps1 First input vertex
		\param ps2 Second input vertex
		\return A non-negative number, the Euclidean distance
	*/
	float getDistance(const HPointSeed &ps1, const HPointSeed &ps2);


	

	/*! Makes the points from \p vertices "visited" if they are situated in the truncated projections of the triangles from the given list \p newTriangles
		\param newTriangles The list of triangles
	*/
	void visitPoints(std::list<HTriangleSeed> &newTriangles);

	/*! Makes the points from \p vertices "visited" if they are situated within the triangle prism of the given \p triangle with the height \p maxSurfaceDepth
		\param triangle The input triangle
	*/
	void visitPoints(HTriangleSeed &triangle);

	/*! Marks \p p as "visited"
		\param p A pointer to a vertex from \p vertices
	*/
	void visitPoint(HPointSeed* p);

	/*! Indicates intersection of two given triangles \p t1 and \p t2 and their nonconformity to one non-selfintersecting surface. 
		\param t1 First input triangle
		\param t2 Second input triangle
		\return True if an approximation of the pyramidal projection of \p t1 intersects the truncated projection of \p t2 in 3D. Otherwise returns false.
	*/
	bool triangles3DIntersection(const common::Triangle<common::Point3<float>> &t1, const common::Triangle<common::Point3<float>> &t2);

	/*! Tests \p triangle3DIntersection() for the given triangle \p t with all triangles from \p triangles
		\param t The input triangle
		\returns True if triangles3DIntersection(t,t1) is true for all t1 from the list \p triangles. Otherwise returns false
		\see triangles3DIntersection
		\see triangles
	*/
	bool triangleMesh3DIntersection(const HTriangleSeed &t);
	
	/*! Checks whether the given triangle \p t is degenerate (at least one of it's angles is too small). Returns true if so, otherwise returns false
		\param t The input triangle
		\return true if the given triangle is degenerate, otherwise returns false
	*/
	bool triangleDegenerate(const HTriangleSeed &t);
	

	/*! Adds the given edge \p e to the list of active edges \p activeEdges. Makes new breakup of \p activeEdges and \p frozenEdges (non-overlapping property) if \p e overlaps any of their members
		\param e The input edge
		\see killOverlappingRegularSegments
		\see activeEdges
	*/
	void addActiveEdge(const HEdgeSeed &e);

	/*! Calculates segment overlapping parameter
		\param ps The input vertex
		\param e The reference edge
		\param t The output parameter value
		\return True and writes into \p t value 0<=\p t<=1 if \p ps can be represented as \p ps=t*\p e.beginPoint+(1-t)*\p e.endPoint. Otherwise returns false
	*/
	bool segmentOverlapParameter(const HPointSeed &ps, const HEdgeSeed &e, float &t);
	
	/*! Deletes all overlapping parts from \p edgeList and \p segmentParts. The number of elements in the both lists can decrease or increase after this procedure
		\param segmentParts First edges list
		\param edgeList Second edges list
	*/
	void killOverlappingRegularSegments(std::list<HEdgeSeed> &segmentParts, std::list<HEdgeSeed> &edgeList);
	
	/*! Alters the coordinates of \p ps in such a way that the angles between the edges of the triangle built from \p e and \p ps, and the elements of \p edgeList are not less than \p maxExcludedAngle. This is achieved by superposition of the nearest edges (sticking them together).
		\param e The input edge
		\param ps The pointer to the input vertex
		\param edgeList The list of the tested edges
		\return True if \p ps was changed, otherwise - false
		\see maxExcludedAngle
	*/
	bool stickToAdjacentEdge(const HEdgeSeed &e, HPointSeed* &ps, std::list<HEdgeSeed> &edgeList);

	/*! Performs \p stickToAdjacentEdge() for \p activeEdges and \p frozenEdges (if the result for \p activeEdges is false), and returns their OR value.
		\param e The input edge
		\param ps The pointer to the input vertex
		\return \p activeEdges() OR \p frozenEdges()
		\see activeEdges
		\see frozenEdges
	*/
	bool excludeSmallAngles(const HEdgeSeed &e, HPointSeed* &ps);
	

	
	
	/*! Calculates the normal vector for the face of the given triangle \p t. Attention: the direction depends on the vertices order
		\param t The input triangle
		\return The normal vector for \p t
	*/
	common::Point3<float> getNormalVector(const common::Triangle<common::Point3<float>> &t);

	/*! Calculates the normal vector for the face of the given triangle \p ts. Attention: the direction depends on the vertices order
		\param ts The input triangle
		\return The normal vector for \p ts
	*/
	common::Point3<float> getNormalVector(const HTriangleSeed &ts);

	/*! Calculates the normal vector of the plane formed by two given vectors \p v1 and \p v2
		\param v1 The first input vector
		\param v2 The second input vector
		\return The normal vector to <v1,v2>
	*/
	common::Point3<float> getNormalVector(const common::Point3<float> v1, const common::Point3<float> v2);
	
	/*! Calculates the element from \p vertices that is nearest to the point of propagation for the given edge \p e. Searches among the visited points if \p checkVisited is true 
		\param e The input vertex
		\param checkVisited A flag. Search among the visited vertices if true
		\return The pointer to the propagated vertex.
	*/
	HPointSeed* getPropagatedVertex(const HEdgeSeed &e, bool checkVisited);
	
	/*! Calculates propagation vectors for the given edges \p e1, \p e2 and \p e3 and modify the edges with them
		\param e1 The first edge
		\param e2 The second edge
		\param e3 The third edge
		\return True if the vectors were successfully calculated and embedded into \p e1,\p e2 and \p e3. Otherwise returns false
	*/
	bool getEdgesPropagations(HEdgeSeed &e1, HEdgeSeed &e2, HEdgeSeed &e3);

	/*! Calculates an approximation of the normal surface vector in point \p p. The surface is defined by the points' cloud within \p vertices. The Procedure is using PCA for the neighborhood of \p p with radius \p windowRadius. The normal vector is defined as the eigenvector with the smallest eigenvalue
		\param p The reference point
		\param windowRadius The radius of the neighborhood
		\return The normal vector
	*/
	common::Point3<float> getSurfaceNormal(common::Point3<float> p, float windowRadius);


	//! Container of the input vertices. Internal realization as a k-DTree
	HVertexContainer *vertices;

	//! List of active edges
	std::list<HEdgeSeed> activeEdges;

	//! List of passive edges
	std::list<HEdgeSeed> frozenEdges;

	//! List of mesh triangles
	std::list<HTriangleSeed> triangles;


	//! The minimal allowed length of the initial triangle's side
	float minInitDistance;

	//! The maximal allowed length of the initial triangle's side
	float maxInitDistance;

	//! The maximal allowed distance between a point and it's node-projection in 3D
	float maxProjectionNodeDistance;

	//! The maximum analyzed depth of the points cloud layer
	float maxSurfaceDepth;

	//! The maximum cos-value of angles between mesh edges that will be excluded during the mesh construction
	float maxExcludedAngle;

	//! The radius of the neighborhood system used for PCA-based normal vector calculation
	float normalNeighborhoodRadius;

	//! The linear proportion between PCA-normal and triangle-based normals
	float faceSurfaceFactor;

	//! The maximum cos-value of the angles between adjacent frozen edges, that allows post-stitching of the edges
	float maxStitchedAngle;

	//! Auxiliary variable. Current polygon square
	float initSquare;

	//! Auxiliary variable. The number of the unvisited vertices
	unsigned int unvisitedCount;
};

} // namespace surfaces


} // namespace methods

#endif //D25_ACTIVE_CONTOURS_HPP_408B8C5F_B876_4B70_AE3C_4B193F9AEED0_