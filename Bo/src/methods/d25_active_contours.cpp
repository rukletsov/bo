
#include "pch.h"
#include <functional>
#include <map>

#include "3rdparty/svd/svd.h"
#include "bo/kdtree.hpp"
#include "bo/blas/blas.hpp"
#include "bo/methods/d25_active_contours.hpp"

namespace bo {
namespace methods {
namespace surfaces {


// HTringleElement implementation.

bool HTriangleElement::operator==(const HTriangleElement &other) const
{
    return (*p1 == *other.p1) && (*p2 == *other.p2) && (*p3 == *other.p3); 
}


// HPointElement implementation.

HPointElement::HPointElement(Vertex v /*= bo::Vector<float,3>(0,0,0)*/)
{
    p = v;
    isVisited = false;
}

bool HPointElement::isNode()
{
    return (adjacentTriangles.size() > 0);
}

float HPointElement::operator[](const size_t t) const
{
    return (t == 0) ? p.x() : ( (t == 1) ? p.y() : p.z() );
}

bool HPointElement::operator==(const HPointElement &other) const
{
    return (p == other.p) && (adjacentTriangles == other.adjacentTriangles) && 
        (isVisited == other.isVisited);
}


// HEdgeElement implementation.

HEdgeElement::HEdgeElement(HPointElement* p1, HPointElement* p2) : p1(p1), p2(p2)
{
}

HEdgeElement::HEdgeElement() : p1(0), p2(0)
{
}

bool HEdgeElement::operator==(const HEdgeElement& other) const
{
    // The propagation direction is not considered.
    return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
}

void HEdgeElement::swap()
{
    HPointElement* tmp = p1;
    p1 = p2;
    p2 = tmp;
}


// Point container implementation.
// Represents a wrapper of HPointElement for utilization in HPointContainer. 

struct HPointContainerItem
{
    // Default constructor.
    HPointContainerItem() : ps(0)
    {
    }

    // Constructor.
    HPointContainerItem(HPointElement* ps) : ps(ps)
    {
    }

    // Access operator.
    inline float operator[](const size_t t) const 
    {
        return (*ps)[t];
    }

    // Comparison operator.
    inline bool operator==(const HPointContainerItem &other) const
    {
        return (*ps) == (*other.ps);
    }

    // Pointer to the reference point element.
    HPointElement* ps;
};


// HPointContainerItem brackets accessor.
inline float point_bac(HPointContainerItem t, size_t k)
{ 
    return t[k]; 
}

// 3D Tree type.
typedef bo::KDTree<3, HPointContainerItem,
    std::pointer_to_binary_function<HPointContainerItem, size_t, float> > D3Tree;

// A general container of vertices.
class HPointContainer
{
public:

    HPointContainer(std::vector<Vector<float,3> >& vertices):
                    tree(D3Tree(std::ptr_fun(point_bac)))
    {
        // Filling in the 3D Tree.
        for (std::vector<Vector<float,3> >::const_iterator itp = vertices.begin();
            itp != vertices.end(); ++itp)
        {
            HPointElement* pe = new HPointElement(*itp);
            
            HPointContainerItem ce(pe);
            tree.insert(ce);

        }
        tree.optimise();
    }

    ~HPointContainer()
    {
    }

    D3Tree tree;
};

// Predicate for closest point with minimal allowed distance.
class PredicateClosestPointBeyondMinDistance
{
public:

    PredicateClosestPointBeyondMinDistance(HPointContainerItem const& searchCenter, float minDistance, 
                                           bool checkNodes, bool checkVisited):
                                           searchCenter(*searchCenter.ps), minDistance(minDistance), 
                                           checkNodes(checkNodes), checkVisited(checkVisited)
    {
    }

    inline bool operator()(HPointContainerItem const& ce) const
    {
        return ((!ce.ps->isVisited) || (checkNodes&&ce.ps->isNode()) || 
                (checkVisited && ce.ps->isVisited)) &&
               ((searchCenter.p-ce.ps->p).eucl_norm() > minDistance);
    }

protected:

    HPointElement searchCenter;

    float minDistance;

    bool checkNodes;

    bool checkVisited;
};

// Predicate for closest point with non-collinearity property.
class PredicateClosestPointNonCollinear
{
public:

    PredicateClosestPointNonCollinear(const HPointContainerItem &ce1, const HPointContainerItem& ce2,
                                      bool checkNodes, bool checkVisited) :
                                      checkNodes(checkNodes), checkVisited(checkVisited),
                                      ce1(ce1), ce2(ce2)
    {
        eps = 0.5;
        this->ce1 = ce1;
        this->ce2 = ce2;
        ba = ce2.ps->p - ce1.ps->p;
        absBa = ba.eucl_norm();
        only_nodes = false;
    }

    inline void check_only_nodes(bool only_nodes)
    {
        this->only_nodes = only_nodes;
    }

    inline bool operator()(HPointContainerItem const& ce) const
    {
        if (only_nodes)
        {
            if (!ce.ps->isNode())
                return false;
        }
        else
        {
            bool isPretender = (!ce.ps->isVisited) || (checkNodes&&ce.ps->isNode()) ||
                               (checkVisited && ce.ps->isVisited);
            if (!isPretender)
                return false;
        }

        if (absBa == 0)
            return false;

        // Collinearity test.
        Vector<float,3> ca = ce.ps->p - ce1.ps->p;;
        double absCa2 = ca.x() * ca.x() + ca.y() * ca.y() + ca.z() * ca.z();

        double scalBaCa = ba * ca;
        double projCaOnBa = scalBaCa / absBa;

        double residual2 = absCa2 - projCaOnBa * projCaOnBa;

        if(residual2 > eps)
            return true;

        else return false;
    }

protected:

    bool checkNodes;

    bool only_nodes;

    bool checkVisited;

    double eps;

    HPointContainerItem ce1;

    HPointContainerItem ce2;

    Vector<float,3> ba;

    double absBa;
};


//General functions


/*! Calculates the normal vector of the plane formed by two given vectors \p v1 and \p v2
    \param v1 The first input vector
    \param v2 The second input vector
    \return The normal vector to <v1,v2>
*/
Vector<float,3> getNormalVector(const Vector<float,3> a, const Vector<float,3> b)
{
    Vector<float,3> p;

    p.x() = a.y() * b.z() - b.y() * a.z();
    p.y() = a.z() * b.x() - b.z() * a.x();
    p.z() = a.x() * b.y() - b.x() * a.y();

    return p;
}

/*! Calculates the normal vector for the face of the given triangle \p t.
    Attention: the direction depends on the vertices order
    \param t The input triangle
    \return The normal vector for \p t
*/
Vector<float, 3> getNormalVector(const Triangle<Vector<float, 3> > &t)
{
    Vector<float,3> a = t.B() - t.A();
    Vector<float,3> b = t.C() - t.A();

    return getNormalVector(a, b);
}

/*! Calculates the normal vector for the face of the given triangle \p ts.
    Attention: the direction depends on the vertices order
    \param ts The input triangle
    \return The normal vector for \p ts
*/
Vector<float,3> getNormalVector(const surfaces::HTriangleElement &ts)
{
    Vector<float,3> a = ts.p2->p - ts.p1->p;
    Vector<float,3> b = ts.p3->p - ts.p1->p;

    return getNormalVector(a, b);
}




// Triangular dipyramid implmentation.


struct TriangularDipyramid
{
    Vector<float,3> vertices[5];
    Triangle<Vector<float,3> > faces[6];
    Vector<float,3> center;
    
    // Triangular dipyramid based on the given triangle with the given base angle.
    static TriangularDipyramid from_triangle_and_angle(const Triangle<Vector<float,3> >& t,
                                                       float baseAngleCos)
    {
        TriangularDipyramid tdp;
        
        // Normal calculation.
        Vector<float,3> z = getNormalVector(t);

        // Sides length.
        float a = float((t.B() - t.C()).eucl_norm());
        float b = float((t.C() - t.A()).eucl_norm());
        float c = float((t.A() - t.B()).eucl_norm());
        
        // Half-perimeter.
        float p = (a + b + c) / 2;

        // Radius of the incircle.
        float r = std::sqrt((p - a) * (p - b) * ( p - c) / p);

        // Center of the incircle.
        tdp.center = (t.A() * a + t.B() * b + t.C() * c) / (a + b + c);

        // Calculate the height of the pyramid such that cos of the angles 
        // between the faces and the base are baseAngleCos.
        float h = r * std::sqrt(1 / (baseAngleCos * baseAngleCos)-1);

        // Height vector.
        z = z / float(z.eucl_norm()) * h;
        
        // Two top-vertices.
        Vector<float,3> p1 = tdp.center + z;
        Vector<float,3> p2 = tdp.center - z;

        // Vertices of the dipyramid.
        tdp.vertices[0] = t.A(); tdp.vertices[1] = t.B(); tdp.vertices[2] = t.C();
        tdp.vertices[3] = p1; tdp.vertices[4] = p2;

        // Faces of the dipyramid.
        tdp.faces[0] = Triangle<Vector<float, 3> >(t.A(), t.B(), p1);
        tdp.faces[1] = Triangle<Vector<float, 3> >(t.B(), t.C(), p1);
        tdp.faces[2] = Triangle<Vector<float, 3> >(t.C(), t.A(), p1);
        tdp.faces[3] = Triangle<Vector<float, 3> >(t.A(), t.B(), p2);
        tdp.faces[4] = Triangle<Vector<float, 3> >(t.B(), t.C(), p2);
        tdp.faces[5] = Triangle<Vector<float, 3> >(t.C(), t.A(), p2);

        return tdp;
    }

    // Triangular dipyramid based on the given triangle with the given height.
    static TriangularDipyramid from_triangle_and_height(const Triangle<Vector<float, 3> >& t,
                                                        float height)
    {
        TriangularDipyramid tdp;

        // Normal calculation.
        Vector<float, 3> z = getNormalVector(t);

        // Sides length.
        float a = float((t.B() - t.C()).eucl_norm());
        float b = float((t.C() - t.A()).eucl_norm());
        float c = float((t.A() - t.B()).eucl_norm());

        // Center of the incircle.
        tdp.center=(t.A() * a + t.B() * b + t.C() * c) / ( a + b + c);

        // Height vector.
        z = z / float(z.eucl_norm()) * height;

        // Two top-vertices.
        Vector<float, 3> p1 = tdp.center + z;
        Vector<float, 3> p2 = tdp.center - z;

        // Vertices of the dipyramid.
        tdp.vertices[0] = t.A(); tdp.vertices[1] = t.B(); tdp.vertices[2] = t.C();
        tdp.vertices[3] = p1; tdp.vertices[4] = p2;

        // Faces of the dipyramid.
        tdp.faces[0] = Triangle<Vector<float, 3> >(t.A(), t.B(), p1);
        tdp.faces[1] = Triangle<Vector<float, 3> >(t.B(), t.C(), p1);
        tdp.faces[2] = Triangle<Vector<float, 3> >(t.C(), t.A(), p1);
        tdp.faces[3] = Triangle<Vector<float, 3> >(t.A(), t.B(), p2);
        tdp.faces[4] = Triangle<Vector<float, 3> >(t.B(), t.C(), p2);
        tdp.faces[5] = Triangle<Vector<float, 3> >(t.C(), t.A(), p2);

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
        for (unsigned i = 0; i < 2; ++i)
        {
            for (unsigned int t = 0; t < 6; ++t)
            {
                bo::Triangle<Vertex > face = tp1.faces[t];
                Vertex norm = getNormalVector(face);

                // Calculating min and max of the projection of the first dipyramid 
                // on the current normal vector.
                float t1MinProj = 0, t1MaxProj = 0;
                tp1.get_min_max_projection(norm, t1MinProj, t1MaxProj);

                // Calculating min and max of the projection of the second dipyramid 
                // on the current normal vector.
                float t2MinProj = 0, t2MaxProj = 0;
                tp2.get_min_max_projection(norm, t2MinProj, t2MaxProj);

                // If the projection intervals do not intersect, 
                // than the convex polygons are separated.
                if(t1MaxProj < t2MinProj + eps || t2MaxProj < t1MinProj + eps)
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

                        Vertex norm = getNormalVector(e1, e2);
                        
                        // For non-zero cross products perform the test.
                        float d = static_cast<float>(norm.eucl_norm());
                        if (d > float_zero_eps)
                        {
                            // Normalize, just for order.
                            norm = norm / d;

                            // Calculating min and max of the projection of the first dipyramid 
                            // on the current normal vector.
                            float t1MinProj = 0, t1MaxProj = 0;
                            tp1.get_min_max_projection(norm, t1MinProj, t1MaxProj);

                            // Calculating min and max of the projection of the second dipyramid 
                            // on the current normal vector.
                            float t2MinProj = 0, t2MaxProj = 0;
                            tp2.get_min_max_projection(norm, t2MinProj, t2MaxProj);

                            // If the projection intervals do not intersect, 
                            // than the convex polygons are separated.
                            if(t1MaxProj < t2MinProj + eps || t2MaxProj < t1MinProj + eps)
                                return false;
                        } 
                    }

        return true;
    }

    // Calculates minimum and maximum of the dipyramid projection on the vector v.
    void get_min_max_projection(const Vector<float,3> &v, float &minProjection, float &maxProjection)
    {
        minProjection = 0;
        maxProjection = 0;

        for (int k = 0; k < 5; ++k)
        {
            float dot = v * vertices[k];
           
            if (k == 0)
            {
                minProjection = maxProjection = dot;
            }
            else
            {
                if (maxProjection < dot)
                    maxProjection = dot;
                if (minProjection > dot)
                    minProjection = dot;
            }
        }
    }
};




//Implementation of D25ActiveContours

D25ActiveContours::D25ActiveContours(float averageFaceSide)
{
    minInitDistance = averageFaceSide;
    
    maxInitDistance = 1.5f * averageFaceSide;
    maxProjectionNodeDistance = 0.8f * averageFaceSide;
    normalNeighborhoodRadius = 0.8f * averageFaceSide;
    maxSurfaceDepth = 0.5f;

    maxExcludedAngle = 0.99f;
    inertialFactor = 0.5;
    tetrahedronBaseAngle = 0.95f;

    vertices = 0;
}

D25ActiveContours::D25ActiveContours(float minInitDistance, float maxInitDistance,
                                     float maxProjectionNodeDistance, float normalNeighborhoodRadius,
                                     float maxSurfaceDepth, float maxExcludedAngle,
                                     float inertialFactor, float tetrahedronBaseAngle):
    minInitDistance(minInitDistance), maxInitDistance(maxInitDistance),
    maxProjectionNodeDistance(maxProjectionNodeDistance), maxSurfaceDepth(maxSurfaceDepth),
    maxExcludedAngle(maxExcludedAngle), normalNeighborhoodRadius(normalNeighborhoodRadius),
    inertialFactor(inertialFactor), tetrahedronBaseAngle(tetrahedronBaseAngle)
{
    vertices = 0;
}

D25ActiveContours::~D25ActiveContours()
{
    delete vertices;
}

inline HPointElement* D25ActiveContours::get_closest_point(const HPointElement &ps,
                                                        bool checkNodes, bool checkVisited)
{
    HPointContainerItem ce(0);

    PredicateClosestPointBeyondMinDistance pred(HPointContainerItem(const_cast<HPointElement*>(&ps)),
                                              minInitDistance,checkNodes,checkVisited);
    std::pair<D3Tree::const_iterator, float> nif = vertices->tree.find_nearest_if(
                HPointContainerItem(const_cast<HPointElement*>(&ps)),maxInitDistance,pred);
    if (nif.first != vertices->tree.end())
        ce=*nif.first;

    return ce.ps;
}

HPointElement* D25ActiveContours::get_closest_min_func_point(const HPointElement &ps1,
    const HPointElement& ps2, bool checkNodes, bool checkVisited)
{
    Vector<float,3> v1 = ps2.p - ps1.p;
    double a = v1.eucl_norm();

    HPointElement mid;
    mid.p = (ps1.p + ps2.p) / 2;

    std::vector<HPointContainerItem> v;
    float const range = maxInitDistance;
    vertices->tree.find_within_range(HPointContainerItem(&mid), range, std::back_inserter(v));

    HPointContainerItem ce(0);
    double func = -1;

    std::vector<HPointContainerItem>::const_iterator it = v.begin();
    while (it != v.end())
    {
        if ((!(*it).ps->isVisited) || (checkNodes&&(*it).ps->isNode()) || 
            (checkVisited&&(*it).ps->isVisited))
        {	
            Vector<float,3> v2 = (*it).ps->p - ps2.p;
            Vector<float,3> v3 = ps1.p - (*it).ps->p;

            double b = v2.eucl_norm();
            double c = v3.eucl_norm();

            if (b < maxInitDistance && c < maxInitDistance && 
                b > minInitDistance && c > minInitDistance)
            {
                double f = std::abs(a-b) + std::abs(a-c);

                if (func == -1 || f < func)
                {
                    func = f;
                    ce = *it;
                }
            }

        } 

        it++;
    }

    return ce.ps;
}


inline HPointElement* D25ActiveContours::get_closest_noncollinear_point(const HPointElement &ps,
    const HPointElement &ps1, const HPointElement &ps2, bool checkNodes, bool checkVisited)
{
    HPointContainerItem ce(0);

    PredicateClosestPointNonCollinear pred(HPointContainerItem(const_cast<HPointElement*>(&ps1)),
                                           HPointContainerItem(const_cast<HPointElement*>(&ps2)),
                                           checkNodes, checkVisited);
    pred.check_only_nodes(true);

    std::pair<D3Tree::const_iterator,float> nif = vertices->tree.find_nearest_if(
        HPointContainerItem(const_cast<HPointElement*>(&ps)), maxProjectionNodeDistance, pred);

    if (nif.first != vertices->tree.end())
    {
        ce = *nif.first;
    }
    else
    {
        pred.check_only_nodes(false);

        nif = vertices->tree.find_nearest_if(HPointContainerItem(const_cast<HPointElement*>(&ps)),
                                             maxProjectionNodeDistance, pred);
        if(nif.first != vertices->tree.end())
        {
            ce = *nif.first;
        }
    }

    return ce.ps;
}

inline float D25ActiveContours::get_distance(const HPointElement &ps1, const HPointElement &ps2)
{
    return (float)(ps1.p - ps2.p).eucl_norm();
}


void D25ActiveContours::model_init()
{
    HPointElement *pps1 = 0, *pps2 = 0, *pps3 = 0;

    D3Tree::mutable_iterator it = vertices->tree.begin();
    while (it != vertices->tree.end())
    {
        // Take the first unvisited point from the point cloud.
        if (!(*it).ps->isVisited)
        {
            pps1 = (*it).ps;

            // Mark it as visited.
            visit_point(pps1);

            // Try to find the closest point that lie out of the sphere
            // with radius minInitialDistance and center in pps1.
            pps2 = get_closest_point(*pps1, true, false);

            // If the point exists:
            if (pps2)
            {
                // Mark it as visited.
                visit_point(pps2);

                // Try to find the point pps3 from the point cloud that makes
                // the triangle [pps1, pps2, pps3] so equilateral as possible.
                pps3 = get_closest_min_func_point(*pps1, *pps2, true, false);

                // If such point exists:
                if (pps3)
                {
                    // Create a new triangle.
                    HTriangleElement tr;
                    tr.p1 = pps1;
                    tr.p2 = pps2;
                    tr.p3 = pps3;

                    // Mark the vertex as visited.
                    visit_point(pps3);

                    // Test if the triangle is well-shaped and does not
                    // intersect the mesh. 
                    if (!triangle_degenerate(tr) && 
                        !triangle_mesh_3d_intersection(tr))
                    {
                        // Create three new edge elements.
                        HEdgeElement e1, e2, e3;
                                                
                        e1.p1 = pps1; e1.p2 = pps2;
                        e2.p1 = pps2; e2.p2 = pps3;
                        e3.p1 = pps3; e3.p2 = pps1;
                        
                        // Try to calculate propagation vectors for them.
                        if (get_edge_propagation(e1, pps3->p) &&
                            get_edge_propagation(e2, pps1->p) &&
                            get_edge_propagation(e3, pps2->p))
                        {
                            // Add the edges to the list of active edges.
                            add_active_edge(e1);
                            add_active_edge(e2);
                            add_active_edge(e3);

                            // Process the points from the point cloud that are in the triangle projection 
                            // (mark as visited). Add the reference of the triangle to its vertices.
                            visit_points(tr);

                            // Add new triangle into the mesh.
                            triangles.push_back(tr);

                            return;
                        }
                    }
                }
            }
        }

        it++;
    }

    return;
}


void D25ActiveContours::model_grow()
{
    if (activeEdges.size() == 0)
        return;

    // Take the first edge from the list.
    HEdgeElement e = *activeEdges.begin();

    // Try to calculate the propagated point from the point cloud.
    HPointElement* pps = get_propagated_vertex(e, false);

    if(pps)
    {
        // Exclude small angles by sticking the propagated
        // point to the ends of adjacent edges.
        exclude_small_angles(e, pps);

        // Create a new triangle.
        HTriangleElement tr;
        tr.p1=e.p1;
        tr.p2=e.p2;
        tr.p3=pps;

        // If the triangle is bad-shaped, put the new edge into
        // the list of passive edges for further processing.
        if(triangle_degenerate(tr))
        {
            frozenEdges.push_back(e);
        }
        // Test the triangle for collision with the mesh. If does not intersect:
        else if(!triangle_mesh_3d_intersection(tr))
        {
            // New edges based on the old one and the propagated point.
            HEdgeElement e1(pps, e.p1);
            HEdgeElement e2(pps, e.p2);
           
            // If the propagation vectors are successfully calculated:
            if(get_edge_propagation(e1, e.p2->p) &&
               get_edge_propagation(e2, e.p1->p))
            {
                // Add two new active edges into the processing list.
                add_active_edge(e1);
                add_active_edge(e2);
               
                // Process the points from the point cloud that are in the triangle projection 
                // (mark as visited). Add the reference of the triangle to its vertices.
                visit_points(tr);	

                // Add new triangle into the mesh.
                triangles.push_back(tr);
            } 
  
        }
        // If the triangle intersects the mesh:
        else 
        {
            // "On-the-fly" stitching.
            edge_stitch(e);
        }
    }
    // If can not calculate the propagated point:
    else 
    {
        frozenEdges.push_back(e);
    }

    // Delete the processed (first) edge from the active edges if it is yet here.
    std::list<HEdgeElement>::iterator it = activeEdges.begin();
    if(it != activeEdges.end() && e == *it)
        activeEdges.erase(it);
}

bool D25ActiveContours::get_edge_propagation(HEdgeElement &e, Vertex origin)
{           
    // The middle point of the edge.
    Vector<float,3> mid = (e.p1->p + e.p2->p) / 2;   

    // The median segment[origin, mid].
    Vector<float,3> median = mid - origin;

    // Propagation norm is sqrt(3)/2 of min distance.
    // Tends to the median length of a equilateral triangle.
    float pnorm = minInitDistance * 0.87f;
    
    // If (inertialFactor >= 1): Inertial edge propagation.
    Vector<float,3> p(0, 0, 0);
    float nmedian = (float)median.eucl_norm();
    if (nmedian != 0)
        p = median / nmedian * pnorm;

    // If (inertialFactor < 1): Tangential PCA-based or adaptive (complex) edge propagation. 
    if (inertialFactor < 1)
    {
        size_t neighbourCount = 0;

        // Surface normal calculation in the middle point of the edge.
        Vector<float,3> midNorm = get_surface_normal(mid, normalNeighborhoodRadius, neighbourCount);

        // Vector parallel to the edge.
        Vector<float,3> v = e.p2->p - e.p1->p;

        // Propagation direction as cross product of v and the normal.
        Vector<float,3> prop = v.cross_product(midNorm);

        // The propagation vector is corrected to be directed outward of the triangle
        // formed by the edge e and the origin point [e, origin].
        // Warning: may not work in regions with extremely high curvature (corners less than 90grad).
        {
            // Check the angle between the median and the propagation vector. 
            float cosa = prop * median;

            // Change the direction of the propagation if it "looks inside".
            if (cosa < 0 ) 
                prop = -prop;
        }

        // Normalize the propagation vector.
        float nprop = (float)prop.eucl_norm();
        if (nprop != 0)
            prop = prop / nprop * pnorm;

        // Tangential PCA-based edge propagation.
        if (inertialFactor >= 0)
        {
            // Linear combination of the inertial and tangential propagations.
            p = (1 - inertialFactor) * prop + inertialFactor * p;
        }
        // Adaptive (complex) edge propagation. Lambda is defined as -inertialFactor.
        else 
        {            
            // Adaptive saturation.
            float alpha = neighbourCount > -inertialFactor ? 1 : neighbourCount / (-inertialFactor);

            // Linear combination of the inertial and tangential propagations.
            p = (1 - alpha) * p + alpha * prop;
        }
    }

    // Insert the calculated propagation into the given edge.
    e.propagationVector = p;

    return true;
}

inline HPointElement* D25ActiveContours::get_propagated_vertex(const HEdgeElement &e,
                                                            bool checkVisited)
{
    HPointElement p;

    // The middle point of the edge.
    Vector<float,3> mid = (e.p1->p + e.p2->p) / 2;

    // Estimate the propagated vertex position.
    p.p = mid + e.propagationVector;

    // Find the closest vertex from the point cloud.
    HPointElement* ps=get_closest_noncollinear_point(p,*e.p1,*e.p2,true, checkVisited);

    return ps;
}

void D25ActiveContours::visit_points( std::list<HTriangleElement> &newTriangles )
{
    std::list<HTriangleElement>::iterator tit = newTriangles.begin();
    while(tit!=newTriangles.end())
    {
        visit_points(*tit);
        ++tit;
    }
}

void D25ActiveContours::visit_points(HTriangleElement &tr)
{
    const float eps = 0.001f;

    // Insert the current triangle in the list of adjacent triangles of 
    // the triangle's nodes.
    tr.p1->adjacentTriangles.push_back(tr);
    tr.p2->adjacentTriangles.push_back(tr);
    tr.p3->adjacentTriangles.push_back(tr);

    // The mass center of the triangle.
    Vector<float,3> mid = (tr.p1->p + tr.p2->p + tr.p3->p) / 3;

    // Define the neighbourhood radius. Find the maximum of maxSurfaceDepth and 
    // the distances from the triangle vertices to its mass centre.
    float tmp, rad = maxSurfaceDepth;
    rad = rad > (tmp = float((tr.p1->p - mid).eucl_norm())) ? rad : tmp;
    rad = rad > (tmp = float((tr.p2->p - mid).eucl_norm())) ? rad : tmp;
    rad = rad > (tmp = float((tr.p3->p - mid).eucl_norm())) ? rad : tmp;

    // Pre-search: choose all points in the range (sphere with radius = rad).
    HPointElement mids;
    mids.p = mid;
    std::vector<HPointContainerItem> v;
    vertices->tree.find_within_range(HPointContainerItem(&mids), rad, std::back_inserter(v));

    // Check if the selected points are inside the triangle prism with 
    // the height = 2*maxSurfaceDepth.
    {
        // Calculate the triangle prism coordinate axes.
        Vector<float,3> X = tr.p2->p - tr.p1->p;
        Vector<float,3> Y = tr.p3->p - tr.p1->p;
        Vector<float,3> Z = getNormalVector(tr); 
        Z = Z / float(Z.eucl_norm()) * maxSurfaceDepth;
        Vector<float,3> O = tr.p1->p;

        // Calculate the triangle prism basis matrix.
        blas::matrix<float> m(3,3);
        m(0,0) = X.x(); m(0,1) = Y.x(); m(0,2) = Z.x();
        m(1,0) = X.y(); m(1,1) = Y.y(); m(1,2) = Z.y();
        m(2,0) = X.z(); m(2,1) = Y.z(); m(2,2) = Z.z();
        
        blas::matrix<float> im(3,3);
        
        // Find the coordinates of the vertices from the range in the calculated
        // coordinate system using the inverse basis matrix. 
        if (blas::invert_matrix(m,im))
        {
            std::vector<HPointContainerItem>::iterator it = v.begin();
            while (it!=v.end())
            {
                if (!it->ps->isVisited)
                {
                    // The coordinates in the old coordinate system.
                    blas::matrix<float> mp(3,1);
                    mp(0,0) = it->ps->p.x() - O.x();
                    mp(1,0) = it->ps->p.y() - O.y();
                    mp(2,0) = it->ps->p.z() - O.z();
                    
                    // Calculate the coordinates in the new (triangle prism)
                    // coordinate system.
                    blas::matrix<float> am = prod(im,mp);

                    float a = am(0,0);
                    float b = am(1,0);
                    float c = am(2,0);

                    // Check if the coordinates are inside the prism.
                    if ((a + b <= 1 + eps) && (a >= -eps) && (b >= -eps) 
                        && (c >= -1) && (c <= 1))
                    {
                        // If so, mark the point as visited.
                        visit_point(it->ps);	
                    }
                }
                ++it;
            }
        }
    }
    
    
    
}

inline void D25ActiveContours::visit_point( HPointElement* p )
{
    if(p->isVisited == false)
    {
        p->isVisited = true;
        --unvisitedCount;
    }
}

inline bool D25ActiveContours::exclude_small_angles( const HEdgeElement &e, HPointElement* &ps )
{
    // Try to connect the given point element ps to the ends of adjacent edges from
    // the lists of active and passive edges (consiquently).
    if (stick_to_adjacent_edge(e, ps, activeEdges) ||
        stick_to_adjacent_edge(e, ps, frozenEdges))
    {
        return true;
    }

    return false;
}

bool D25ActiveContours::stick_to_adjacent_edge(const HEdgeElement &e, HPointElement* &ps,
                                               std::list<HEdgeElement> &edgeList)
{  
    std::list<HEdgeElement>::const_iterator it = edgeList.begin();
    while (it != edgeList.end())
    {
        HEdgeElement ee = *it;

        // Skip the base edge e.
        if (((e.p1->p == ee.p1->p) && (e.p2->p == ee.p2->p)) ||
            ((e.p1->p == ee.p2->p) && (e.p2->p == ee.p1->p)))
        {
            ++it;
            continue;
        }

        // If at least one of the new edges [e.p1, ps] or [e.p2, ps] is coinciding with
        // the current edge, finish the procedure (the vertex ps is already sticked).
        if (((e.p1->p == ee.p1->p || e.p2->p == ee.p1->p) && (ps->p == ee.p2->p)) ||
            ((e.p1->p == ee.p2->p || e.p2->p == ee.p2->p) && (ps->p == ee.p1->p)))
        {
            return true;
        }

        // Sticking procedure.
        {
            // The common vertex of e and ee.
            Vector<float,3> p;

            // Define the adjacent edge.
            if (((p = e.p1->p) == ee.p1->p) || (e.p1->p == ee.p2->p) ||
                ((p = e.p2->p) == ee.p1->p) || (e.p2->p == ee.p2->p))
            {
                // Order the verices of the adgacent edge.
                if (p == ee.p2->p)
                    ee.swap();
                
                // Vector [ps, p]. 
                Vector<float,3> v1 = ps->p - p;
                double normV1 = v1.eucl_norm();

                // Vector [ee.p2, p].
                Vector<float,3> v2 = ee.p2->p - p;
                double normV2 = v2.eucl_norm();

                if(normV1==0 || normV2==0)
                    return false;

                // Calculate cos of the angle between the adjacent vectors.
                double cosa = v1 * v2 / float(normV1 * normV2);

                // If the angle is less than minimal allowed, perform sticking.
                if(cosa > maxExcludedAngle)
                {
                    ps = ee.p2;
                    return true;
                }

            }
        }
        

        ++it;
    }

    return false;
}


Mesh D25ActiveContours::build_mesh(std::vector<Vertex> &v)
{
    // Initialize the point cloud.
    set_vertices(v);

    // Build triangles while the growing step is possible.
    while (grow_step())
    { }

    // Generate and return the mesh.
    return get_mesh();
}

Mesh D25ActiveContours::build_mesh()
{
    // Clean and initialize auxiliary structures.
    prepare();

    // Build triangles while the growing step is possible.
    while (grow_step())
    { }

    // Generate and return the mesh.
    return get_mesh();
}

inline bool D25ActiveContours::triangle_mesh_3d_intersection(const HTriangleElement &t)
{
    // Calculate the mass center of the triangle.
    Triangle<Vector<float,3> > t1(t.p1->p, t.p2->p, t.p3->p);
    HPointElement mass((t1.A() + t1.B() + t1.C())/3);

    // Calculate the neighbourhood of the mass center.
    std::vector<HPointContainerItem> neighbours;
    float const range = 3*minInitDistance;
    vertices->tree.find_within_range(HPointContainerItem(&mass), range, std::back_inserter(neighbours));

    // For all nodes from the neighbourhood check their adjacent triangles for 
    // intersection with the given triangle.
    std::vector<HPointContainerItem>::const_iterator it = neighbours.begin();
    while(it != neighbours.end())
    {
        std::list<HTriangleElement>::const_iterator tit = it->ps->adjacentTriangles.begin();
        while(tit != it->ps->adjacentTriangles.end())
        {
            Triangle<Vector<float, 3> > t2(tit->p1->p ,tit->p2->p, tit->p3->p);

            if (triangles_3d_intersection(t1, t2))
                return true;

            ++tit;
        }
        ++it;
    }

    return false;
}

void D25ActiveContours::edge_stitch(HEdgeElement e )
{
    // Try to calculate the propagation point from the point cloud.
    HPointElement* pps1 = get_propagated_vertex(e, true);

    bool isStitched = false;

    // If the propagated point exists:
    if (pps1)
    {
        // Create a new triangle based on the given edge and the propagated point.
        Triangle<Vector<float, 3> > t1(e.p1->p, e.p2->p, pps1->p);

        // Find all adjacent edges to the given one in the list of passive edges.
        for(std::list<HEdgeElement>::iterator ite = frozenEdges.begin(); ite != frozenEdges.end(); ++ite)
        {
            HEdgeElement ee = *ite;

            if (e == ee) continue;

            bool b11 = false, b12 = false, b21 = false, b22 = false;

            // Test for adjacency.
            // Suppress warning C4706 using the comparison with true.
            if ( ((b11 = (e.p1 == ee.p1)) == true) || ((b12 = (e.p1 == ee.p2)) == true) ||
                 ((b21 = (e.p2 == ee.p1)) == true) || ((b22 = (e.p2 == ee.p2)) == true) )
            {
                // If the edge is adjacent try to calculate its propagation point
                // from the point cloud.
                HPointElement* pps2 = get_propagated_vertex(ee, true);
                
                // If the such point exists:
                if (pps2)
                {
                    // Create a new triangle based on this adjacent edge and the
                    // propagated point
                    Triangle<Vector<float, 3> > t2(ee.p1->p, ee.p2->p, pps2->p);

                    // If these two triangles are concurrent, create one stitch triangle based
                    // on the considered adjacent edges.
                    if (triangles_3d_intersection(t1, t2))
                    {
                        // Changing everything to b11 condition to order the directions
                        // of the edges.
                        if (b12)
                        {
                            ee.swap();
                        }
                        else if (b21)
                        {
                            e.swap();
                        }
                        else if (b22)
                        {
                            e.swap();
                            ee.swap();
                        }

                        // Vertices of the new triangle.
                        HPointElement *p1, *p2, *p3;
                        p1 = e.p1;
                        p2 = e.p2;
                        p3 = ee.p2;

                        // Edge stitching: new edge of the stitch triangle.
                        HEdgeElement ne(p2, p3);                  

                        // Calculate the propagation vector for the stitch edge.
                        if (get_edge_propagation(ne, p1->p))
                        {
                            // Create the stitch triangle.
                            HTriangleElement tr;
                            tr.p1 = p1;
                            tr.p2 = p2;
                            tr.p3 = p3;
 
                            // Test for collisions of the triangle with the mesh.
                            if (!triangle_mesh_3d_intersection(tr))
                            {
                                // Delete the current frozen neighboring edge.
                                frozenEdges.erase(ite);

                                // Add a new active edge.
                                add_active_edge(ne);

                                // Add new triangle to the mesh.
                                triangles.push_back(tr);

                                isStitched=true;

                                break;
                            } 
                        } 
                    } 
                } // if(pps2)
            }
        }
    } // if(pps1)

    // Add the edge to the frozen edges, if it wasn't stitched on this step.
    if (!isStitched)
        frozenEdges.push_back(e);
}

Vertex D25ActiveContours::get_surface_normal(Vertex p, float windowRadius, std::size_t &neighbourCount)
{
    // Select points from the neighborhood.
    HPointElement ps;
    ps.p = p;
    std::vector<HPointContainerItem> neighbours;
    float const range = windowRadius;
    vertices->tree.find_within_range(HPointContainerItem(&ps), range,
                                     std::back_inserter(neighbours));

    size_t pointCount = neighbours.size();
    neighbourCount = pointCount;

    // Important.
    if (pointCount == 0)
        return Vertex(0, 0, 0);

    // Perform the Principal Component Analysis (PCA):
    
    Vector<float,3> mean(0,0,0);

    // The mean calculation.
    std::vector<HPointContainerItem>::const_iterator itp = neighbours.begin();
    while(itp != neighbours.end())
    {
        Vector<float,3> pp = (*itp).ps->p;
        mean = mean + pp;

        ++itp;
    }
    mean = mean / (float)pointCount;
   
    // The extended covariance matrix preparation [COV | Ex].
    double** A = new double*[6]; 
    for (int i = 0; i < 6; ++i)
    {
        A[i] = new double[3];
        A[i][0] = A[i][1] = A[i][2] = 0;
    }

    double* S2 = new double[3];

    // The covariance matrix calculation, A[0,0 - 3,3].
    {
        itp = neighbours.begin();
        while (itp != neighbours.end())
        {
            Vector<float,3> pp = (*itp).ps->p;

            A[0][0] += (pp.x()-mean.x()) * (pp.x()-mean.x()); 
            A[1][0] += (pp.y()-mean.y()) * (pp.x()-mean.x()); 
            A[2][0] += (pp.z()-mean.z()) * (pp.x()-mean.x()); 
            
            A[0][1] += (pp.x()-mean.x()) * (pp.y()-mean.y()); 
            A[1][1] += (pp.y()-mean.y()) * (pp.y()-mean.y()); 
            A[2][1] += (pp.z()-mean.z()) * (pp.y()-mean.y());
            
            A[0][2] += (pp.x()-mean.x()) * (pp.z()-mean.z()); 
            A[1][2] += (pp.y()-mean.y()) * (pp.z()-mean.z()); 
            A[2][2] += (pp.z()-mean.z()) * (pp.z()-mean.z());

            ++itp;
        }

        for (int i = 0; i < 3; ++i)
        {
            A[i][0] /= pointCount;
            A[i][1] /= pointCount;
            A[i][2] /= pointCount;
        }
    }
 
    // Vector dimension.
    int n = 3;

    // Calculate eigenvectors using SVD decomposition.
    svd(A,S2,n);

    // The index of the minimal element.
    int minIndex = (S2[0] < S2[1]) ? (S2[0] < S2[2] ? 0 : 2) : (S2[1] < S2[2] ? 1 : 2);

    // The eigenvector the corresponds to the minimal eigenvalue. This vector is taken
    // as an approximation the point cloud surface normal.
    Vector<float,3> v((float)A[3][minIndex], (float)A[4][minIndex], (float)A[5][minIndex]);

    // Release the resources.
    for (int i = 0; i < 6; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] S2;

    // Return the normalized normal vector/
    return v / float(v.eucl_norm());
}



std::vector<HPointElement> D25ActiveContours::get_vertices()
{   
    std::vector<HPointElement> verts;

    // Compose a vector from the tree-based container.
    D3Tree::mutable_iterator it = vertices->tree.begin();
    while (it != vertices->tree.end())
    {
        HPointContainerItem c = *it;
        verts.push_back(*c.ps);
        ++it;
    }

    return verts;
}

const std::list<HEdgeElement>* D25ActiveContours::get_active_edges()
{
    return &activeEdges;
}

const std::list<HEdgeElement>* D25ActiveContours::get_frozen_edges()
{
    return &frozenEdges;
}

const std::list<HTriangleElement>* D25ActiveContours::get_triangles()
{
    return &triangles;
}



bool D25ActiveContours::triangles_3d_intersection(const Triangle<Vertex> &t1,
                                                  const Triangle<Vertex> &t2)
{
    // Define the minimal error value.
    const float eps = 0.001f;
    // Bound the angle from zero by eps.
    float alpha = tetrahedronBaseAngle < eps ? eps : tetrahedronBaseAngle;
    
    // Create two dypiramids from the given triangles and the base angle alpha.
    TriangularDipyramid tdp1 = TriangularDipyramid::from_triangle_and_angle(t1, alpha);
    TriangularDipyramid tdp2 = TriangularDipyramid::from_triangle_and_angle(t2, alpha);

    // Test them for intersection.
    return tdp1.intersects(tdp2);
   
}

bool D25ActiveContours::triangle_degenerate(const HTriangleElement &t)
{
    // Adjacent vectors of the triangle angles.
    // Angle a:
    Vector<float,3> v1a = t.p2->p - t.p1->p;
    Vector<float,3> v1b = t.p3->p - t.p1->p;
    // Angle b:
    Vector<float,3> v2a = t.p1->p - t.p2->p;
    Vector<float,3> v2b = t.p3->p - t.p2->p; 
    // Angle c:
    Vector<float,3> v3a = t.p1->p - t.p3->p;
    Vector<float,3> v3b = t.p2->p - t.p3->p;
 
    // Calculate cos of the angles a, b and c.
    float cos1 = v1a * v1b / float((v1a.eucl_norm() * v1b.eucl_norm()));
    float cos2 = v2a * v2b / float((v2a.eucl_norm() * v2b.eucl_norm()));
    float cos3 = v3a * v3b / float((v3a.eucl_norm() * v3b.eucl_norm()));

    // Compare the angles with the minimal allowed value. 
    if ((std::fabs(cos1) > maxExcludedAngle) ||
        (std::fabs(cos2) > maxExcludedAngle) ||
        (std::fabs(cos3) > maxExcludedAngle))
        return true;
    else return false;
}

void D25ActiveContours::set_vertices( std::vector<Vertex> &v )
{
    // Create and fill in a vertex container.
    if (vertices)
        delete vertices;
    vertices=new HPointContainer(v);

    // Clean and initialize auxiliary structures.
    prepare();
}

bool D25ActiveContours::grow_step()
{	
    // If all the points have been processed and the list of the active edges is empty:
    if (unvisitedCount == 0 && activeEdges.size() == 0)
    {
        // Perform the post-stitch step and check if new active or passive edges were
        // added during it (if so, return true).
        std::size_t frozenBefore = frozenEdges.size();
        post_stitch();
        return (activeEdges.size() > 0) || (frozenEdges.size() < frozenBefore);
    }
    // If there are active edges in the list:
    else if (activeEdges.size() > 0)
    {
        // Perform a new growing step.
        model_grow();
    }
    // If there are unprocessed points and no active edges in the list
    // perform the initialization step. 
    else model_init();

    return true;
}

void D25ActiveContours::post_stitch()
{
    // TODO: complete this step if needed.
    // Purpose: it performs mutual stitch of the frozenEdges based on some rule R(edge1, edge2). 
}


Mesh D25ActiveContours::get_mesh()
{
    //Construct a new mesh.
    Mesh m(triangles.size());

    // Create a reference map (from the local triangles nodes to the mesh vertices).
    std::map<HPointElement*,size_t> mymap;
    
    // Fill in the mesh.
    std::list<HTriangleElement>::const_iterator itt = triangles.begin();
    while (itt != triangles.end())
    {
        for (int j = 0; j < 3; ++j)
        {
            // Take one of three triangle's vertices. 
            HPointElement* ps = (j == 0) ? (itt->p1) : (j == 1 ? itt->p2 : itt->p3);

            // If the triangle node is not yet in the map:
            if (mymap.find(ps) == mymap.end())
            {	
                // Add the vertex/node into the mesh and into the reference map.
                size_t ind = m.add_vertex(Mesh::Vertex(ps->p.x(), ps->p.y(), ps->p.z()));
                mymap[ps] = ind;
            }
        }

        // Add the face from the processed triangle nodes.
        size_t A = mymap[itt->p1];
        size_t B = mymap[itt->p2];
        size_t C = mymap[itt->p3];

        m.add_face(Mesh::Face(A,B,C));

        ++itt;
    }

    return m;
}

void D25ActiveContours::prepare()
{
    // Cleaning all the auxiliary containers.
    activeEdges.clear();
    frozenEdges.clear();
    triangles.clear();

    // Initialize the counter.
    unvisitedCount = static_cast<unsigned>(vertices->tree.size());
}

void D25ActiveContours::add_active_edge( HEdgeElement &e )
{
    bool isUnique = true;

    std::list<HEdgeElement>::iterator it = activeEdges.begin();
    
    // Check if the given edge is unique. TODO: optimize!
    while(it != activeEdges.end())
    {
        if (e == *it)
        {
            isUnique = false;
            
            // Two equal active edges kill themselves. 
            activeEdges.erase(it);
            
            break;
        }

        ++it;
    }

    if (isUnique)
        activeEdges.push_back(e);
}

} // namespace surfaces
} // namespace methods
} // namespace bo
