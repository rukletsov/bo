
#include "bo/pch.h"
#include <functional>
#include <map>

#include "3rdparty/svd/svd.h"
#include "bo/core/kdtree.hpp"
#include "bo/math/blas_extensions.hpp"
#include "bo/surfaces/d25_active_contours.hpp"

namespace bo {
namespace surfaces {

namespace detail {

// Point container implementation.
// Represents a wrapper of PointElement for utilization in PointContainer.
struct PointContainerItem
{
    // Default constructor.
    PointContainerItem() : ps(0)
    { }

    // Constructor.
    PointContainerItem(bo::surfaces::detail::PointElement* ps) : ps(ps)
    { }

    // Access operator.
    inline float operator[](const size_t t) const
    {
        return (*ps)[t];
    }

    // Comparison operator.
    inline bool operator==(const PointContainerItem &other) const
    {
        return (*ps) == (*other.ps);
    }

    // Pointer to the reference point element.
    bo::surfaces::detail::PointElement* ps;
};


// HPointContainerItem brackets accessor.
inline float point_bac(PointContainerItem t, size_t k)
{ 
    return t[k]; 
}

// 3D Tree type.
typedef bo::KDTree<3, PointContainerItem,
    std::pointer_to_binary_function<PointContainerItem, size_t, float> > D3Tree;

// A general container of vertices.
class PointContainer
{
public:
    PointContainer(std::vector<bo::surfaces::Vertex> &vertices):
                   tree(D3Tree(std::ptr_fun(point_bac)))
    {
        // Filling in the 3D Tree.
        for (std::vector<bo::surfaces::Vertex>::const_iterator itp = vertices.begin();
            itp != vertices.end(); ++itp)
        {
            bo::surfaces::detail::PointElement* pe = new bo::surfaces::detail::PointElement(*itp);
            
            PointContainerItem ce(pe);
            tree.insert(ce);

        }
        tree.optimise();
    }

    ~PointContainer()
    { }

    D3Tree tree;
};

// Predicate for closest point with minimal allowed distance.
class PredicateClosestPointBeyondMinDistance
{
public:

    PredicateClosestPointBeyondMinDistance(PointContainerItem const& searchCenter,
            float minDistance, bool checkNodes, bool checkVisited):
        searchCenter(*searchCenter.ps), minDistance(minDistance),
        checkNodes(checkNodes), checkVisited(checkVisited)
    { }

    inline bool operator()(PointContainerItem const& ce) const
    {
        return ((!ce.ps->is_visited) || (checkNodes&&ce.ps->is_node()) ||
                (checkVisited && ce.ps->is_visited)) &&
               ((searchCenter.p-ce.ps->p).euclidean_norm_d() > minDistance);
    }

protected:
    bo::surfaces::detail::PointElement searchCenter;
    float minDistance;
    bool checkNodes;
    bool checkVisited;
};

// Predicate for closest point with non-collinearity property.
class PredicateClosestPointNonCollinear
{
public:
    PredicateClosestPointNonCollinear(const PointContainerItem &ce1,
            const PointContainerItem& ce2, bool checkNodes, bool checkVisited) :
        checkNodes(checkNodes), checkVisited(checkVisited), ce1(ce1), ce2(ce2)
    {
        eps = 0.5;
        this->ce1 = ce1;
        this->ce2 = ce2;
        ba = ce2.ps->p - ce1.ps->p;
        absBa = ba.euclidean_norm_d();
        only_nodes = false;
    }

    inline void check_only_nodes(bool only_nodes)
    {
        this->only_nodes = only_nodes;
    }

    inline bool operator()(PointContainerItem const& ce) const
    {
        if (only_nodes)
        {
            if (!ce.ps->is_node())
                return false;
        }
        else
        {
            bool isPretender = (!ce.ps->is_visited) || (checkNodes&&ce.ps->is_node()) ||
                               (checkVisited && ce.ps->is_visited);
            if (!isPretender)
                return false;
        }

        if (absBa == 0)
            return false;

        // Collinearity test.
        bo::surfaces::Vertex ca = ce.ps->p - ce1.ps->p;;
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

    PointContainerItem ce1;
    PointContainerItem ce2;
    bo::surfaces::Vertex ba;
    double absBa;
};

} // namespace detail


//Implementation of D25ActiveContours
D25ActiveContours::D25ActiveContours(float average_face_side)
{
    min_init_distance_ = average_face_side;
    
    max_init_distance_ = 1.5f * average_face_side;
    max_projection_node_distance_ = 0.8f * average_face_side;
    normal_neighborhood_radius_ = 0.8f * average_face_side;
    max_surface_depth_ = 0.5f;

    max_excluded_angle_ = 0.99f;
    inertial_factor_ = 0.5;
    tetrahedron_base_angle_ = 0.95f;

    vertices_ = 0;
}

D25ActiveContours::D25ActiveContours(float min_init_distance, float max_init_distance,
                                     float max_projection_node_distance,
                                     float normal_neighborhood_radius,
                                     float max_surface_depth, float max_excluded_angle,
                                     float inertial_factor, float tetrahedron_base_angle):
    min_init_distance_(min_init_distance), max_init_distance_(max_init_distance),
    max_projection_node_distance_(max_projection_node_distance), max_surface_depth_(max_surface_depth),
    max_excluded_angle_(max_excluded_angle), normal_neighborhood_radius_(normal_neighborhood_radius),
    inertial_factor_(inertial_factor), tetrahedron_base_angle_(tetrahedron_base_angle)
{
    vertices_ = 0;
}

D25ActiveContours::~D25ActiveContours()
{
    delete vertices_;
}

detail::PointElement* D25ActiveContours::get_closest_point(
        const detail::PointElement &ps, bool checkNodes, bool checkVisited)
{
    detail::PointContainerItem ce(0);

    detail::PredicateClosestPointBeyondMinDistance pred(detail::PointContainerItem(
            const_cast<detail::PointElement*>(&ps)), min_init_distance_, checkNodes, checkVisited);
    std::pair<detail::D3Tree::const_iterator, float> nif = vertices_->tree.find_nearest_if(
            detail::PointContainerItem(const_cast<detail::PointElement*>(&ps)),max_init_distance_,pred);
    if (nif.first != vertices_->tree.end())
        ce = *nif.first;

    return ce.ps;
}

detail::PointElement* D25ActiveContours::get_closest_min_func_point(const detail::PointElement &ps1,
        const detail::PointElement &ps2, bool checkNodes, bool checkVisited)
{
    Vertex v1 = ps2.p - ps1.p;
    double a = v1.euclidean_norm_d();

    detail::PointElement mid;
    mid.p = (ps1.p + ps2.p) / 2;

    std::vector<detail::PointContainerItem> v;
    float const range = max_init_distance_;
    vertices_->tree.find_within_range(detail::PointContainerItem(&mid), range, std::back_inserter(v));

    detail::PointContainerItem ce(0);
    double func = -1;

    std::vector<detail::PointContainerItem>::const_iterator it = v.begin();
    while (it != v.end())
    {
        if ((!(*it).ps->is_visited) || (checkNodes&&(*it).ps->is_node()) ||
            (checkVisited && (*it).ps->is_visited))
        {	
            Vertex v2 = (*it).ps->p - ps2.p;
            Vertex v3 = ps1.p - (*it).ps->p;

            double b = v2.euclidean_norm_d();
            double c = v3.euclidean_norm_d();

            if (b < max_init_distance_ && c < max_init_distance_ &&
                b > min_init_distance_ && c > min_init_distance_)
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

detail::PointElement* D25ActiveContours::get_closest_noncollinear_point(
        const detail::PointElement &ps, const detail::PointElement &ps1,
        const detail::PointElement &ps2, bool checkNodes, bool checkVisited)
{
    detail::PointContainerItem ce(0);

    detail::PredicateClosestPointNonCollinear pred(detail::PointContainerItem(
        const_cast<detail::PointElement*>(&ps1)),
        detail::PointContainerItem(const_cast<detail::PointElement*>(&ps2)), checkNodes, checkVisited);
    pred.check_only_nodes(true);

    std::pair<detail::D3Tree::const_iterator,float> nif = vertices_->tree.find_nearest_if(
        detail::PointContainerItem(const_cast<detail::PointElement*>(&ps)),
        max_projection_node_distance_, pred);

    if (nif.first != vertices_->tree.end())
    {
        ce = *nif.first;
    }
    else
    {
        pred.check_only_nodes(false);

        nif = vertices_->tree.find_nearest_if(detail::PointContainerItem(
            const_cast<detail::PointElement*>(&ps)), max_projection_node_distance_, pred);
        if(nif.first != vertices_->tree.end())
        {
            ce = *nif.first;
        }
    }

    return ce.ps;
}

inline
float D25ActiveContours::get_distance(const detail::PointElement &ps1,
                                      const detail::PointElement &ps2)
{
    return (float)(ps1.p - ps2.p).euclidean_norm_d();
}

void D25ActiveContours::model_init()
{
    detail::PointElement *pps1 = 0, *pps2 = 0, *pps3 = 0;

    detail::D3Tree::mutable_iterator it = vertices_->tree.begin();
    while (it != vertices_->tree.end())
    {
        // Take the first unvisited point from the point cloud.
        if (!(*it).ps->is_visited)
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
                    detail::TriangleElement tr;
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
                        detail::EdgeElement e1, e2, e3;
                                                
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
                            triangles_.push_back(tr);

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
    if (active_edges_.size() == 0)
        return;

    // Take the first edge from the list.
    detail::EdgeElement e = *active_edges_.begin();

    // Try to calculate the propagated point from the point cloud.
    detail::PointElement* pps = get_propagated_vertex(e, false);

    if(pps)
    {
        // Exclude small angles by sticking the propagated
        // point to the ends of adjacent edges.
        exclude_small_angles(e, pps);

        // Create a new triangle.
        detail::TriangleElement tr;
        tr.p1=e.p1;
        tr.p2=e.p2;
        tr.p3=pps;

        // If the triangle is bad-shaped, put the new edge into
        // the list of passive edges for further processing.
        if(triangle_degenerate(tr))
        {
            frozen_edges_.push_back(e);
        }
        // Test the triangle for collision with the mesh. If does not intersect:
        else if(!triangle_mesh_3d_intersection(tr))
        {
            // New edges based on the old one and the propagated point.
            detail::EdgeElement e1(pps, e.p1);
            detail::EdgeElement e2(pps, e.p2);
           
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
                triangles_.push_back(tr);
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
        frozen_edges_.push_back(e);
    }

    // Delete the processed (first) edge from the active edges if it is yet here.
    std::list<detail::EdgeElement>::iterator it = active_edges_.begin();
    if(it != active_edges_.end() && e == *it)
        active_edges_.erase(it);
}

bool D25ActiveContours::get_edge_propagation(detail::EdgeElement &e, Vertex origin)
{
    // The middle point of the edge.
    Vertex mid = (e.p1->p + e.p2->p) / 2;

    // The median segment[origin, mid].
    Vertex median = mid - origin;

    // Propagation norm is sqrt(3)/2 of min distance.
    // Tends to the median length of a equilateral triangle.
    float pnorm = min_init_distance_ * 0.87f;
    
    // If (inertialFactor >= 1): Inertial edge propagation.
    Vertex p(0, 0, 0);
    float nmedian = (float)median.euclidean_norm_d();
    if (nmedian != 0)
        p = median / nmedian * pnorm;

    // If (inertialFactor < 1): Tangential PCA-based or adaptive (complex) edge propagation. 
    if (inertial_factor_ < 1)
    {
        size_t neighbourCount = 0;

        // Surface normal calculation in the middle point of the edge.
        Vertex midNorm = get_surface_normal(mid, normal_neighborhood_radius_, neighbourCount);

        // Vector parallel to the edge.
        Vertex v = e.p2->p - e.p1->p;

        // Propagation direction as cross product of v and the normal.
        Vertex prop = v.cross_product(midNorm);

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
        float nprop = (float)prop.euclidean_norm_d();
        if (nprop != 0)
            prop = prop / nprop * pnorm;

        // Tangential PCA-based edge propagation.
        if (inertial_factor_ >= 0)
        {
            // Linear combination of the inertial and tangential propagations.
            p = (1 - inertial_factor_) * prop + inertial_factor_ * p;
        }
        // Adaptive (complex) edge propagation. Lambda is defined as -inertialFactor.
        else 
        {            
            // Adaptive saturation.
            float alpha = neighbourCount > -inertial_factor_ ? 1 : neighbourCount / (-inertial_factor_);

            // Linear combination of the inertial and tangential propagations.
            p = (1 - alpha) * p + alpha * prop;
        }
    }

    // Insert the calculated propagation into the given edge.
    e.propagation_vector = p;

    return true;
}

detail::PointElement* D25ActiveContours::get_propagated_vertex(const detail::EdgeElement &e,
                                                               bool checkVisited)
{
    detail::PointElement p;

    // The middle point of the edge.
    Vertex mid = (e.p1->p + e.p2->p) / 2;

    // Estimate the propagated vertex position.
    p.p = mid + e.propagation_vector;

    // Find the closest vertex from the point cloud.
    detail::PointElement* ps = get_closest_noncollinear_point(p, *e.p1, *e.p2, true,
                                                              checkVisited);

    return ps;
}

void D25ActiveContours::visit_points(std::list<detail::TriangleElement>& newTriangles)
{
    std::list<detail::TriangleElement>::iterator tit = newTriangles.begin();
    while(tit != newTriangles.end())
    {
        visit_points(*tit);
        ++tit;
    }
}

void D25ActiveContours::visit_points(detail::TriangleElement &tr)
{
    const float eps = 0.001f;

    // Insert the current triangle in the list of adjacent triangles of 
    // the triangle's nodes.
    tr.p1->adjacent_triangles.push_back(tr);
    tr.p2->adjacent_triangles.push_back(tr);
    tr.p3->adjacent_triangles.push_back(tr);

    // The mass center of the triangle.
    Vertex mid = (tr.p1->p + tr.p2->p + tr.p3->p) / 3;

    // Define the neighbourhood radius. Find the maximum of max_surface_depth and
    // the distances from the triangle vertices to its mass centre.
    float tmp, rad = max_surface_depth_;
    rad = rad > (tmp = float((tr.p1->p - mid).euclidean_norm_d())) ? rad : tmp;
    rad = rad > (tmp = float((tr.p2->p - mid).euclidean_norm_d())) ? rad : tmp;
    rad = rad > (tmp = float((tr.p3->p - mid).euclidean_norm_d())) ? rad : tmp;

    // Pre-search: choose all points in the range (sphere with radius = rad).
    detail::PointElement mids;
    mids.p = mid;
    std::vector<detail::PointContainerItem> v;
    vertices_->tree.find_within_range(detail::PointContainerItem(&mids), rad, std::back_inserter(v));

    // Check if the selected points are inside the triangle prism with 
    // the height = 2*maxSurfaceDepth.
    {
        // Calculate the triangle prism coordinate axes.
        Vertex X = tr.p2->p - tr.p1->p;
        Vertex Y = tr.p3->p - tr.p1->p;
        Vertex Z = tr.normal_vector();
        Z = Z / float(Z.euclidean_norm_d()) * max_surface_depth_;
        Vertex O = tr.p1->p;

        // Calculate the triangle prism basis matrix.
        math::matrix<float> m(3,3);
        m(0,0) = X.x(); m(0,1) = Y.x(); m(0,2) = Z.x();
        m(1,0) = X.y(); m(1,1) = Y.y(); m(1,2) = Z.y();
        m(2,0) = X.z(); m(2,1) = Y.z(); m(2,2) = Z.z();
        
        math::matrix<float> im(3,3);
        
        // Find the coordinates of the vertices from the range in the calculated
        // coordinate system using the inverse basis matrix. 
        if (math::invert_matrix(m,im))
        {
            std::vector<detail::PointContainerItem>::iterator it = v.begin();
            while (it!=v.end())
            {
                if (!it->ps->is_visited)
                {
                    // The coordinates in the old coordinate system.
                    math::matrix<float> mp(3,1);
                    mp(0,0) = it->ps->p.x() - O.x();
                    mp(1,0) = it->ps->p.y() - O.y();
                    mp(2,0) = it->ps->p.z() - O.z();
                    
                    // Calculate the coordinates in the new (triangle prism)
                    // coordinate system.
                    math::matrix<float> am = prod(im,mp);

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

inline
void D25ActiveContours::visit_point(detail::PointElement* p)
{
    if(p->is_visited == false)
    {
        p->is_visited = true;
        --unvisited_count_;
    }
}

inline
bool D25ActiveContours::exclude_small_angles(const detail::EdgeElement &e, detail::PointElement* &ps)
{
    // Try to connect the given point element ps to the ends of adjacent edges from
    // the lists of active and passive edges (consiquently).
    if (stick_to_adjacent_edge(e, ps, active_edges_) ||
        stick_to_adjacent_edge(e, ps, frozen_edges_))
    {
        return true;
    }

    return false;
}

bool D25ActiveContours::stick_to_adjacent_edge(const detail::EdgeElement &e,
        detail::PointElement* &ps, std::list<detail::EdgeElement> &edgeList)
{  
    std::list<detail::EdgeElement>::const_iterator it = edgeList.begin();
    while (it != edgeList.end())
    {
        detail::EdgeElement ee = *it;

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
            Vertex p;

            // Define the adjacent edge.
            if (((p = e.p1->p) == ee.p1->p) || (e.p1->p == ee.p2->p) ||
                ((p = e.p2->p) == ee.p1->p) || (e.p2->p == ee.p2->p))
            {
                // Order the verices of the adgacent edge.
                if (p == ee.p2->p)
                    ee.swap();
                
                // Vector [ps, p]. 
                Vertex v1 = ps->p - p;
                double normV1 = v1.euclidean_norm_d();

                // Vector [ee.p2, p].
                Vertex v2 = ee.p2->p - p;
                double normV2 = v2.euclidean_norm_d();

                if(normV1==0 || normV2==0)
                    return false;

                // Calculate cos of the angle between the adjacent vectors.
                double cosa = v1 * v2 / float(normV1 * normV2);

                // If the angle is less than minimal allowed, perform sticking.
                if(cosa > max_excluded_angle_)
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

void D25ActiveContours::prepare()
{
    // Cleaning all the auxiliary containers.
    active_edges_.clear();
    frozen_edges_.clear();
    triangles_.clear();

    // Initialize the counter.
    unvisited_count_ = static_cast<unsigned>(vertices_->tree.size());
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

bool D25ActiveContours::grow_step()
{
    // If all the points have been processed and the list of the active edges is empty:
    if (unvisited_count_ == 0 && active_edges_.size() == 0)
    {
        return false;
    }
    // If there are active edges in the list:
    else if (active_edges_.size() > 0)
    {
        // Perform a new growing step.
        model_grow();
    }
    // If there are unprocessed points and no active edges in the list
    // perform the initialization step.
    else model_init();

    return true;
}

inline bool D25ActiveContours::triangle_mesh_3d_intersection(const detail::TriangleElement &t)
{
    // Calculate the mass center of the triangle.
    Triangle<Vertex> t1(t.p1->p, t.p2->p, t.p3->p);
    detail::PointElement mass((t1.A() + t1.B() + t1.C())/3);

    // Calculate the neighbourhood of the mass center.
    std::vector<detail::PointContainerItem> neighbours;
    float const range = 3 * max_init_distance_;
    vertices_->tree.find_within_range(detail::PointContainerItem(&mass), range,
                                      std::back_inserter(neighbours));

    // For all nodes from the neighbourhood check their adjacent triangles for 
    // intersection with the given triangle.
    std::vector<detail::PointContainerItem>::const_iterator it = neighbours.begin();
    while(it != neighbours.end())
    {
        std::list<detail::TriangleElement>::const_iterator tit = it->ps->adjacent_triangles.begin();
        while(tit != it->ps->adjacent_triangles.end())
        {
            Triangle<Vertex> t2(tit->p1->p ,tit->p2->p, tit->p3->p);

            if (triangles_3d_intersection(t1, t2))
                return true;

            ++tit;
        }
        ++it;
    }

    return false;
}

void D25ActiveContours::edge_stitch(detail::EdgeElement e)
{
    // Try to calculate the propagation point from the point cloud.
    detail::PointElement* pps1 = get_propagated_vertex(e, true);

    bool isStitched = false;

    // If the propagated point exists:
    if (pps1)
    {
        // Create a new triangle based on the given edge and the propagated point.
        Triangle<Vertex> t1(e.p1->p, e.p2->p, pps1->p);

        // Find all adjacent edges to the given one in the list of passive edges.
        for(std::list<detail::EdgeElement>::iterator ite = frozen_edges_.begin();
            ite != frozen_edges_.end(); ++ite)
        {
            detail::EdgeElement ee = *ite;

            if (e == ee) continue;

            bool b11 = false, b12 = false, b21 = false, b22 = false;

            // Test for adjacency.
            // Suppress warning C4706 using the comparison with true.
            if ( ((b11 = (e.p1 == ee.p1)) == true) || ((b12 = (e.p1 == ee.p2)) == true) ||
                 ((b21 = (e.p2 == ee.p1)) == true) || ((b22 = (e.p2 == ee.p2)) == true) )
            {
                // If the edge is adjacent try to calculate its propagation point
                // from the point cloud.
                detail::PointElement* pps2 = get_propagated_vertex(ee, true);
                
                // If the such point exists:
                if (pps2)
                {
                    // Create a new triangle based on this adjacent edge and the
                    // propagated point
                    Triangle<Vertex> t2(ee.p1->p, ee.p2->p, pps2->p);

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
                        detail::PointElement *p1, *p2, *p3;
                        p1 = e.p1;
                        p2 = e.p2;
                        p3 = ee.p2;

                        // Edge stitching: new edge of the stitch triangle.
                        detail::EdgeElement ne(p2, p3);

                        // Calculate the propagation vector for the stitch edge.
                        if (get_edge_propagation(ne, p1->p))
                        {
                            // Create the stitch triangle.
                            detail::TriangleElement tr;
                            tr.p1 = p1;
                            tr.p2 = p2;
                            tr.p3 = p3;
 
                            // Test for collisions of the triangle with the mesh.
                            if (!triangle_mesh_3d_intersection(tr))
                            {
                                // Delete the current frozen neighboring edge.
                                frozen_edges_.erase(ite);

                                // Add a new active edge.
                                add_active_edge(ne);

                                // Add new triangle to the mesh.
                                triangles_.push_back(tr);

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
        frozen_edges_.push_back(e);
}

Vertex D25ActiveContours::get_surface_normal(Vertex p, float windowRadius, std::size_t &neighbourCount)
{
    // Select points from the neighborhood.
    detail::PointElement ps;
    ps.p = p;
    std::vector<detail::PointContainerItem> neighbours;
    float const range = windowRadius;
    vertices_->tree.find_within_range(detail::PointContainerItem(&ps), range,
                                      std::back_inserter(neighbours));

    size_t pointCount = neighbours.size();
    neighbourCount = pointCount;

    // Important.
    if (pointCount == 0)
        return Vertex(0, 0, 0);

    // Perform the Principal Component Analysis (PCA):
    
    Vertex mean(0,0,0);

    // The mean calculation.
    std::vector<detail::PointContainerItem>::const_iterator itp = neighbours.begin();
    while(itp != neighbours.end())
    {
        Vertex pp = (*itp).ps->p;
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
            Vertex pp = (*itp).ps->p;

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
    Vertex v((float)A[3][minIndex], (float)A[4][minIndex], (float)A[5][minIndex]);

    // Release the resources.
    for (int i = 0; i < 6; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] S2;

    // Return the normalized normal vector/
    return v / float(v.euclidean_norm_d());
}

std::vector<detail::PointElement> D25ActiveContours::get_vertices()
{   
    std::vector<detail::PointElement> verts;

    // Compose a vector from the tree-based container.
    detail::D3Tree::mutable_iterator it = vertices_->tree.begin();
    while (it != vertices_->tree.end())
    {
        detail::PointContainerItem c = *it;
        verts.push_back(*c.ps);
        ++it;
    }

    return verts;
}

const std::list<detail::EdgeElement>* D25ActiveContours::get_active_edges()
{
    return &active_edges_;
}

const std::list<detail::EdgeElement>* D25ActiveContours::get_frozen_edges()
{
    return &frozen_edges_;
}

const std::list<detail::TriangleElement>* D25ActiveContours::get_triangles()
{
    return &triangles_;
}

bool D25ActiveContours::triangles_3d_intersection(const Triangle<Vertex> &t1,
                                                  const Triangle<Vertex> &t2)
{
    // Define the minimal error value.
    const float eps = 0.001f;
    // Bound the angle from zero by eps.
    float alpha = tetrahedron_base_angle_ < eps ? eps : tetrahedron_base_angle_;
    
    // Create two dypiramids from the given triangles and the base angle alpha.
    detail::TriangularDipyramid tdp1 = detail::TriangularDipyramid::from_triangle_and_angle(t1, alpha);
    detail::TriangularDipyramid tdp2 = detail::TriangularDipyramid::from_triangle_and_angle(t2, alpha);

    // Test them for intersection.
    return tdp1.intersects(tdp2);
   
}

bool D25ActiveContours::triangle_degenerate(const detail::TriangleElement &t)
{
    // Adjacent vectors of the triangle angles.
    // Angle a:
    Vertex v1a = t.p2->p - t.p1->p;
    Vertex v1b = t.p3->p - t.p1->p;
    // Angle b:
    Vertex v2a = t.p1->p - t.p2->p;
    Vertex v2b = t.p3->p - t.p2->p;
    // Angle c:
    Vertex v3a = t.p1->p - t.p3->p;
    Vertex v3b = t.p2->p - t.p3->p;
 
    // Calculate cos of the angles a, b and c.
    float cos1 = v1a * v1b / float((v1a.euclidean_norm_d() * v1b.euclidean_norm_d()));
    float cos2 = v2a * v2b / float((v2a.euclidean_norm_d() * v2b.euclidean_norm_d()));
    float cos3 = v3a * v3b / float((v3a.euclidean_norm_d() * v3b.euclidean_norm_d()));

    // Compare the angles with the minimal allowed value. 
    if ((std::fabs(cos1) > max_excluded_angle_) ||
        (std::fabs(cos2) > max_excluded_angle_) ||
        (std::fabs(cos3) > max_excluded_angle_))
        return true;
    else return false;
}

void D25ActiveContours::set_vertices(std::vector<Vertex> &v)
{
    // Create and fill in a vertex container.
    if (vertices_)
        delete vertices_;
    vertices_ = new detail::PointContainer(v);

    // Clean and initialize auxiliary structures.
    prepare();
}

Mesh D25ActiveContours::get_mesh()
{
    //Construct a new mesh.
    Mesh m(triangles_.size());

    // Create a reference map (from the local triangles nodes to the mesh vertices).
    std::map<detail::PointElement*,size_t> mymap;
    
    // Fill in the mesh.
    std::list<detail::TriangleElement>::const_iterator itt = triangles_.begin();
    while (itt != triangles_.end())
    {
        for (int j = 0; j < 3; ++j)
        {
            // Take one of three triangle's vertices. 
            detail::PointElement* ps = (j == 0) ? (itt->p1) : (j == 1 ? itt->p2 : itt->p3);

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

void D25ActiveContours::add_active_edge(detail::EdgeElement &e)
{
    bool isUnique = true;

    std::list<detail::EdgeElement>::iterator it = active_edges_.begin();
    
    // Check if the given edge is unique. TODO: optimize!
    while(it != active_edges_.end())
    {
        if (e == *it)
        {
            isUnique = false;
            
            // Two equal active edges kill themselves. 
            active_edges_.erase(it);
            
            break;
        }

        ++it;
    }

    if (isUnique)
        active_edges_.push_back(e);
}

} // namespace surfaces
} // namespace bo
