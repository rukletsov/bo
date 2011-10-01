
#include "pch.h"
#include <functional>
#include <map>

#include "svd/svd.h"
#include "kdtree++/kdtree.hpp"
#include "common/blas/blas.hpp"
#include "common/methods/d25_active_contours.hpp"

//Uncomment for equal-square triangles propagation  
//#define USE_EQUAL_SQUARE_PROPAGATION

using namespace common;

namespace common {
namespace methods {

namespace {

/*! \struct HTreeElement
\brief HPointSeed wrapper for HVertexContainer utilization 
*/
struct HContainerElement
{
    /*! Default constructor */
    HContainerElement():ps(0){}

    /*! Constructor */
    HContainerElement(surfaces::HPointSeed* ps):ps(ps){}

    /*! Access operator */
    inline float operator [] (const size_t t) const 
    {
        return (*ps)[t];
    }

    /*! Comparison operator */
    inline bool operator == (const HContainerElement &other) const
    {
        return (*ps)==(*other.ps);
    }

    /*! Pointer to an elementary vertex item*/
    surfaces::HPointSeed* ps;
};


//Brackets accessor
inline float bac( HContainerElement t, size_t k ) { return t[k]; }

//3D Tree type
typedef KDTree::KDTree<3, HContainerElement,
    std::pointer_to_binary_function<HContainerElement,size_t,float> > D3Tree;

//3D Tree Wrapper
class HVertexContainer
{
public:

    HVertexContainer(std::vector<Vector<float,3> >& vertices):
        tree(D3Tree(std::ptr_fun(bac)))
    {		
        linear.resize(vertices.size());

        //Filling in the 3D Tree
        int cnt=0;
        for(std::vector<Vector<float,3> >::const_iterator itp = vertices.begin();
            itp != vertices.end(); ++itp)
        {
            linear[cnt].isNode=false;
            linear[cnt].isVisited=false;
            linear[cnt].p=*itp;

            HContainerElement ce(&linear[cnt]);
            tree.insert(ce);

            ++cnt;
        }
        tree.optimise();
    }

    ~HVertexContainer()
    {
    }

    D3Tree tree;
    std::vector<surfaces::HPointSeed> linear;
};

//Predicate for closest point with minimal allowed distance
class PredicateClosestPointWithMinDistance
{
public:
    PredicateClosestPointWithMinDistance(HContainerElement const& searchCenter, float minDistance, bool checkNodes, bool checkVisited):searchCenter(*searchCenter.ps),minDistance(minDistance),checkNodes(checkNodes),checkVisited(checkVisited)
    {
    }
    inline bool operator()( HContainerElement const& ce ) const
    {
        return ((!ce.ps->isVisited)||(checkNodes&&ce.ps->isNode)||(checkVisited&&ce.ps->isVisited))&&((searchCenter.p-ce.ps->p).eucl_norm()>minDistance);
    }
protected:
    surfaces::HPointSeed searchCenter;
    float minDistance;
    bool checkNodes;
    bool checkVisited;
};

//Predicate for closest point with non-collinearity property
class PredicateClosestPointNonCollinear
{
public:
    PredicateClosestPointNonCollinear(const HContainerElement &ce1, const HContainerElement& ce2, bool checkNodes, bool checkVisited):ce1(ce1),ce2(ce2),checkNodes(checkNodes),checkVisited(checkVisited)
    {
        eps=0.5;
        this->ce1=ce1;
        this->ce2=ce2;
        ba=ce2.ps->p-ce1.ps->p;
        absBa=ba.eucl_norm();
        only_nodes=false;
    }
    inline void check_only_nodes(bool only_nodes)
    {
        this->only_nodes=only_nodes;
    }
    inline bool operator()( HContainerElement const& ce ) const
    {
        if(only_nodes)
        {
            if(!ce.ps->isNode)return false;
        }
        else
        {
            bool isPretender=(!ce.ps->isVisited)||(checkNodes&&ce.ps->isNode)||(checkVisited&&ce.ps->isVisited);
            if(!isPretender)return false;
        }

        if(absBa==0)return false;

        //Collinearity test
        Vector<float,3> ca=ce.ps->p-ce1.ps->p;;
        double absCa2=ca.x()*ca.x()+ca.y()*ca.y()+ca.z()*ca.z();

        double scalBaCa=ba*ca;
        double projCaOnBa=scalBaCa/absBa;

        double residual2=absCa2-projCaOnBa*projCaOnBa;

        if(residual2>eps)return true;
        else return false;
    }
protected:
    bool checkNodes;
    bool only_nodes;
    bool checkVisited;
    double eps;
    HContainerElement ce1;
    HContainerElement ce2;
    Vector<float,3> ba;
    double absBa;
};

/*! Calculates the normal vector of the plane formed by two given vectors \p v1 and \p v2
    \param v1 The first input vector
    \param v2 The second input vector
    \return The normal vector to <v1,v2>
*/
Vector<float,3> getNormalVector( const Vector<float,3> a, const Vector<float,3> b )
{
    Vector<float,3> p;

    p.x()=a.y()*b.z()-b.y()*a.z();
    p.y()=a.z()*b.x()-b.z()*a.x();
    p.z()=a.x()*b.y()-b.x()*a.y();

    return p;
}

/*! Calculates the normal vector for the face of the given triangle \p t. Attention: the direction depends on the vertices order
    \param t The input triangle
    \return The normal vector for \p t
*/
Vector<float, 3> getNormalVector(const Triangle<Vector<float, 3> > &t)
{
    Vector<float,3> a=t.B()-t.A();
    Vector<float,3> b=t.C()-t.A();

    return getNormalVector(a,b);
}

/*! Calculates the normal vector for the face of the given triangle \p ts. Attention: the direction depends on the vertices order
    \param ts The input triangle
    \return The normal vector for \p ts
*/
Vector<float,3> getNormalVector(const surfaces::HTriangleSeed &ts)
{
    Vector<float,3> a=ts.p2->p-ts.p1->p;
    Vector<float,3> b=ts.p3->p-ts.p1->p;

    return getNormalVector(a,b);
}

//Triangular dipyramid
struct TriangularDipyramid
{
    Vector<float,3> vertices[5];
    Triangle<Vector<float,3> > faces[6];
    Vector<float,3> center;
    
    //Triangular dipyramid based on the given triangle with the given base angle
    static TriangularDipyramid from_triangle_and_angle(const Triangle<Vector<float,3> >& t,
                                                       float baseAngleCos)
    {
        TriangularDipyramid tdp;
        
        //Normal calculation.
        Vector<float,3> z = getNormalVector(t);

        //Sides length.
        float a=float((t.B()-t.C()).eucl_norm());
        float b=float((t.C()-t.A()).eucl_norm());
        float c=float((t.A()-t.B()).eucl_norm());
        
        //Half-perimeter.
        float p=(a+b+c)/2;

        //Radius of the incircle.
        float r=sqrt((p-a)*(p-b)*(p-c)/p);

        //Center of the incircle.
        tdp.center=(t.A()*a+t.B()*b+t.C()*c)/(a+b+c);

        //Height of the pyramid such that cos of the angles between the faces
        //and the base are baseAngleCos.
        float h=r*sqrt(1/(baseAngleCos*baseAngleCos)-1);

        //Height vector
        z=z/float(z.eucl_norm())*h;
        
        //Two top-vertices
        Vector<float,3> p1=tdp.center+z;
        Vector<float,3> p2=tdp.center-z;

        //Vertices of the dipyramid
        tdp.vertices[0]=t.A(); tdp.vertices[1]=t.B(); tdp.vertices[2]=t.C();
        tdp.vertices[3]=p1; tdp.vertices[4]=p2;

        //Faces of the dipyramid
        tdp.faces[0]=Triangle<Vector<float, 3> >(t.A(), t.B(), p1);
        tdp.faces[1]=Triangle<Vector<float, 3> >(t.B(), t.C(), p1);
        tdp.faces[2]=Triangle<Vector<float, 3> >(t.C(), t.A(), p1);
        tdp.faces[3]=Triangle<Vector<float, 3> >(t.A(), t.B(), p2);
        tdp.faces[4]=Triangle<Vector<float, 3> >(t.B(), t.C(), p2);
        tdp.faces[5]=Triangle<Vector<float, 3> >(t.C(), t.A(), p2);

        return tdp;
    }

    //Triangular dipyramid based on the given triangle with the given height
    static TriangularDipyramid from_triangle_and_height(const Triangle<Vector<float, 3> >& t,
                                                        float height)
    {
        TriangularDipyramid tdp;

        //Normal calculation.
        Vector<float,3> z=getNormalVector(t);

        //Sides length.
        float a=float((t.B()-t.C()).eucl_norm());
        float b=float((t.C()-t.A()).eucl_norm());
        float c=float((t.A()-t.B()).eucl_norm());

        //Center of the incircle.
        tdp.center=(t.A()*a+t.B()*b+t.C()*c)/(a+b+c);

        //Height vector
        z=z/float(z.eucl_norm())*height;

        //Two top-vertices
        Vector<float,3> p1=tdp.center+z;
        Vector<float,3> p2=tdp.center-z;

        //Vertices of the dipyramid
        tdp.vertices[0]=t.A(); tdp.vertices[1]=t.B(); tdp.vertices[2]=t.C();
        tdp.vertices[3]=p1; tdp.vertices[4]=p2;

        //Faces of the dipyramid
        tdp.faces[0]=Triangle<Vector<float,3> >(t.A(), t.B(), p1);
        tdp.faces[1]=Triangle<Vector<float,3> >(t.B(), t.C(), p1);
        tdp.faces[2]=Triangle<Vector<float,3> >(t.C(), t.A(), p1);
        tdp.faces[3]=Triangle<Vector<float,3> >(t.A(), t.B(), p2);
        tdp.faces[4]=Triangle<Vector<float,3> >(t.B(), t.C(), p2);
        tdp.faces[5]=Triangle<Vector<float,3> >(t.C(), t.A(), p2);

        return tdp;
    }

    //Check intersection with another triangular dipyramid.
    //Based on the Separatiing Plane Theorem.
    bool intersects(TriangularDipyramid& other)
    {
        float eps=0.01f;

        TriangularDipyramid tp1=*this, tp2=other;

        for(unsigned i=0; i<2; ++i)
        {
            for(unsigned int t=0; t<6; ++t)
            {
                Triangle<Vector<float,3> > face=tp1.faces[t];
                Vector<float,3> norm=getNormalVector(face);
                
                //Calculating min and max of the projection of the first dipyramid 
                //on the current normal vector
                float t1MinProj=0, t1MaxProj=0;
                for(int k=0; k<5; ++k)
                {
                    Vector<float,3> v=tp1.vertices[k];
                    float dot=norm*v;
                    if(k==0)
                    {
                        t1MaxProj=t1MinProj=dot;
                    }
                    else
                    {
                        if(t1MaxProj<dot)t1MaxProj=dot;
                        if(t1MinProj>dot)t1MinProj=dot;
                    }
                }

                //Calculating min and max of the projection of the second dipyramid 
                //on the current normal vector
                float t2MinProj=0, t2MaxProj=0;
                for(int k=0; k<5; ++k)
                {
                    Vector<float,3> v=tp2.vertices[k];
                    float dot=norm*v;
                    if(k==0)
                    {
                        t2MaxProj=t2MinProj=dot;
                    }
                    else
                    {
                        if(t2MaxProj<dot)t2MaxProj=dot;
                        if(t2MinProj>dot)t2MinProj=dot;
                    }
                }

                //If the projection intervals do not intersect, than the convex polygons are separated
                if(t1MaxProj<t2MinProj+eps||t2MaxProj<t1MinProj+eps)return false;
            }

            //swap the dipyramids
            TriangularDipyramid tmp=tp1;
            tp1=tp2;
            tp2=tmp;
        }

        return true;
    }
};

} //anonymous namespace


namespace surfaces {

D25ActiveContours::D25ActiveContours(float minFaceInitSize)
{
    minInitDistance=minFaceInitSize;
    
    maxInitDistance=1.5f*minInitDistance;
    maxProjectionNodeDistance=0.5f*minInitDistance;
    normalNeighborhoodRadius=maxInitDistance;
    maxSurfaceDepth=0.5f;

    maxExcludedAngle=0.90f;
    maxStitchedAngle=-0.90f;
    faceSurfaceFactor=0.0;

    tetrahedronBaseAngle=0.7f;

    vertices=0;
}

D25ActiveContours::D25ActiveContours(float minInitDistance, float maxInitDistance, float maxProjectionNodeDistance,
                                     float normalNeighborhoodRadius, float maxSurfaceDepth, float maxExcludedAngle,
                                     float maxStitchedAngle, float faceSurfaceFactor, float tetrahedronBaseAngle):
                                     minInitDistance(minInitDistance), 
                                     maxInitDistance(maxInitDistance),
                                     maxProjectionNodeDistance(maxProjectionNodeDistance), 
                                     normalNeighborhoodRadius(normalNeighborhoodRadius),
                                     maxSurfaceDepth(maxSurfaceDepth),
                                     maxExcludedAngle(maxExcludedAngle),
                                     maxStitchedAngle(maxStitchedAngle),
                                     faceSurfaceFactor(faceSurfaceFactor),
                                     tetrahedronBaseAngle(tetrahedronBaseAngle)
{
    vertices=0;
}

D25ActiveContours::~D25ActiveContours()
{
    delete vertices;
}

inline HPointSeed* D25ActiveContours::get_closest_point( const HPointSeed &ps, bool checkNodes, bool checkVisited )
{
    HContainerElement ce(0);

    PredicateClosestPointWithMinDistance pred(HContainerElement(const_cast<HPointSeed*>(&ps)),minInitDistance,checkNodes,checkVisited);
    std::pair<D3Tree::const_iterator,float> nif = vertices->tree.find_nearest_if(HContainerElement(const_cast<HPointSeed*>(&ps)),maxInitDistance,pred);
    if(nif.first!=vertices->tree.end())ce=*nif.first;

    return ce.ps;
}

HPointSeed* D25ActiveContours::get_closest_min_func_point( const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited )
{
    Vector<float,3> v1=ps2.p-ps1.p;
    double a=v1.eucl_norm();

    HPointSeed mid;
    mid.p=(ps1.p+ps2.p)/2;

    //Pre-search: choose all points in the range
    std::vector<HContainerElement> v;
    float const range = maxInitDistance;
    vertices->tree.find_within_range(HContainerElement(&mid), range, std::back_inserter(v));

    HContainerElement ce(0);
    double func=-1;

    std::vector<HContainerElement>::const_iterator it = v.begin();
    while(it!=v.end())
    {
        if((!(*it).ps->isVisited)||(checkNodes&&(*it).ps->isNode)||(checkVisited&&(*it).ps->isVisited))
        {	
            Vector<float,3> v2=(*it).ps->p-ps2.p;
            Vector<float,3> v3=ps1.p-(*it).ps->p;

            double b=v2.eucl_norm();
            double c=v3.eucl_norm();

            if(b<maxInitDistance&&c<maxInitDistance&&b>minInitDistance&&c>minInitDistance)
            {
                double f=abs(a-b)+abs(a-c);

                if(func==-1||f<func)
                {
                    func=f;
                    ce=*it;
                }
            }

        } 

        it++;
    }

    return ce.ps;
}


inline HPointSeed* D25ActiveContours::get_closest_noncollinear_point(const HPointSeed &ps, const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited)
{
    HContainerElement ce(0);

    PredicateClosestPointNonCollinear pred(HContainerElement(const_cast<HPointSeed*>(&ps1)),HContainerElement(const_cast<HPointSeed*>(&ps2)),checkNodes,checkVisited);
    pred.check_only_nodes(true);
    std::pair<D3Tree::const_iterator,float> nif = vertices->tree.find_nearest_if(HContainerElement(const_cast<HPointSeed*>(&ps)),maxProjectionNodeDistance,pred);
    if(nif.first!=vertices->tree.end())ce=*nif.first;
    else
    {
        pred.check_only_nodes(false);
        nif = vertices->tree.find_nearest_if(HContainerElement(const_cast<HPointSeed*>(&ps)),maxProjectionNodeDistance,pred);
        if(nif.first!=vertices->tree.end())ce=*nif.first;
    }

    return ce.ps;
}




inline float D25ActiveContours::get_distance( const HPointSeed &ps1, const HPointSeed &ps2 )
{
    return (float)(ps1.p-ps2.p).eucl_norm();
}


void D25ActiveContours::model_init()
{
    HPointSeed *pps1=0,*pps2=0,*pps3=0;

    D3Tree::mutable_iterator it = vertices->tree.begin();
    while(it!=vertices->tree.end())
    {
        if(!(*it).ps->isVisited)
        {
            pps1=(*it).ps;

            visit_point(pps1);

            pps2=get_closest_point(*pps1, true, false);

            if(pps2)
            {
                visit_point(pps2);

                pps3=get_closest_min_func_point(*pps1, *pps2, true, false);

                if(pps3)
                {
                    HTriangleSeed tr;
                    tr.p1=pps1;
                    tr.p2=pps2;
                    tr.p3=pps3;

                    visit_point(pps3);

                    if(!triangle_degenerate(tr)&&!triangle_mesh_3d_intersection(tr))
                    {
                        HEdgeSeed e1, e2, e3;

                        e1.p1=e3.p2=pps1;
                        e2.p1=e1.p2=pps2;
                        e3.p1=e2.p2=pps3;

                        if(get_edges_propagations(e1,e2,e3))
                        {
                            add_active_edge(e1);
                            add_active_edge(e2);
                            add_active_edge(e3);

                            visit_points(tr);

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
    if(activeEdges.size()==0)return;

    HEdgeSeed e=*activeEdges.begin();

    HPointSeed* pps=get_propagated_vertex(e, false);

    if(pps)
    {
        //Exclude small angles
        exclude_small_angles(e, pps);

        HTriangleSeed tr;
        tr.p1=e.p1;
        tr.p2=e.p2;
        tr.p3=pps;

        if(triangle_degenerate(tr))
        {
            frozenEdges.push_back(e);
        }
        //Test for triangles collisions
        else if(!triangle_mesh_3d_intersection(tr))
        {
            //New edges based on an old one and the propagated point
            HEdgeSeed e1(e.p2,pps);
            HEdgeSeed e2(pps,e.p1);
            HEdgeSeed e3=e; //fictive

            //If the propagations are successfully calculated
            if(get_edges_propagations(e1,e2,e3))
            {
                //Add two new active edges
                add_active_edge(e1);
                add_active_edge(e2);

                //Mark visited and under projection points
                visit_points(tr);	

                //Add new triangle to the mesh
                triangles.push_back(tr);
            } 
  
        }
        else 
        {
            //"On-the-fly" stitching
            edge_stitch(e);
        }
    }
    else 
    {
        frozenEdges.push_back(e);
    }

    //Delete the processed edge from the active edges if it is yet here
    std::list<HEdgeSeed>::iterator it=activeEdges.begin();
    if(it!=activeEdges.end()&&e==*it)
        activeEdges.erase(it);
}

bool D25ActiveContours::get_edges_propagations( HEdgeSeed &e1, HEdgeSeed &e2, HEdgeSeed &e3 )
{
    if(e1.p1!=e3.p2 || e2.p1!=e1.p2 || e3.p1!=e2.p2)return false;

    Vector<float,3> mid=(e1.p1->p+e2.p1->p+e3.p1->p)/3;

    Vector<float,3> v1=e1.p2->p-e1.p1->p;
    Vector<float,3> v2=e2.p2->p-e2.p1->p;
    Vector<float,3> v3=e3.p2->p-e3.p1->p;

    Vector<float,3> propagationE1,propagationE2,propagationE3;

    //PCA-based approximation 
    if(faceSurfaceFactor!=0)
    {
        //PCA-based surface normal calculation
        Vector<float,3> midNorm=get_surface_normal(mid,normalNeighborhoodRadius);

        //Propagation directions. Surface normal calculation
        propagationE1=getNormalVector(v1,midNorm);
        propagationE2=getNormalVector(v2,midNorm);
        propagationE3=getNormalVector(v3,midNorm);

        //Outer directions correction
        Vector<float,3> medianE1=(e1.p1->p+e1.p2->p)/2-mid;
        Vector<float,3> medianE2=(e2.p1->p+e2.p2->p)/2-mid;
        Vector<float,3> medianE3=(e3.p1->p+e3.p2->p)/2-mid;

        float cosMP1=(medianE1*propagationE1)/float(medianE1.eucl_norm()*propagationE1.eucl_norm());
        float cosMP2=(medianE2*propagationE2)/float(medianE2.eucl_norm()*propagationE2.eucl_norm());
        float cosMP3=(medianE3*propagationE3)/float(medianE3.eucl_norm()*propagationE3.eucl_norm());

        if(cosMP1<0)propagationE1=propagationE1*(-1);
        if(cosMP2<0)propagationE2=propagationE2*(-1);
        if(cosMP3<0)propagationE3=propagationE3*(-1);

        //Mixture with the PCA-face propagations
        propagationE1=(propagationE1/float(propagationE1.eucl_norm()))*faceSurfaceFactor+(medianE1/float(medianE1.eucl_norm()))*(1-faceSurfaceFactor);
        propagationE2=(propagationE2/float(propagationE2.eucl_norm()))*faceSurfaceFactor+(medianE2/float(medianE2.eucl_norm()))*(1-faceSurfaceFactor);
        propagationE3=(propagationE3/float(propagationE3.eucl_norm()))*faceSurfaceFactor+(medianE3/float(medianE3.eucl_norm()))*(1-faceSurfaceFactor);
    }
    //Triangle face based approximation
    else
    {
        propagationE1=(e1.p1->p+e1.p2->p)/2-mid;
        propagationE2=(e2.p1->p+e2.p2->p)/2-mid;
        propagationE3=(e3.p1->p+e3.p2->p)/2-mid;
    }


    //Normalization
    propagationE1=propagationE1/float(propagationE1.eucl_norm());
    propagationE2=propagationE2/float(propagationE2.eucl_norm());
    propagationE3=propagationE3/float(propagationE3.eucl_norm());

#ifdef USE_EQUAL_SQUARE_PROPAGATION

    //3D propagation triangle square calculation
    float a=v1.getNorm();
    float b=v2.getNorm();
    float c=v3.getNorm();
    float pp=(a+b+c)/2;
    //Heron formula
    initSquare=sqrt(pp*(pp-a)*(pp-b)*(pp-c));

    //New propagations lengths
    float cosE1=(propagationE1*v1)/a;
    float cosE2=(propagationE2*v2)/b;
    float cosE3=(propagationE3*v3)/c;
    float normPE1=2*initSquare/(a*sqrt(1-cosE1*cosE1));
    float normPE2=2*initSquare/(b*sqrt(1-cosE2*cosE2));
    float normPE3=2*initSquare/(c*sqrt(1-cosE3*cosE3));

#else

    float normPE1 = (maxInitDistance+minInitDistance)/2;
    float normPE2 = (maxInitDistance+minInitDistance)/2;
    float normPE3 = (maxInitDistance+minInitDistance)/2;

#endif

    propagationE1=propagationE1*normPE1;
    propagationE2=propagationE2*normPE2;
    propagationE3=propagationE3*normPE3;

    e1.propagationVector=propagationE1;
    e2.propagationVector=propagationE2;
    e3.propagationVector=propagationE3;

    return true;
}


inline HPointSeed* D25ActiveContours::get_propagated_vertex(const HEdgeSeed &e, bool checkVisited)
{
    HPointSeed p;

    Vector<float,3> mid=(e.p1->p+e.p2->p)/2;

    p.p=mid+e.propagationVector;

    HPointSeed* ps=get_closest_noncollinear_point(p,*e.p1,*e.p2,true, checkVisited);

    return ps;
}


void D25ActiveContours::visit_points( std::list<HTriangleSeed> &newTriangles )
{
    std::list<HTriangleSeed>::iterator tit = newTriangles.begin();
    while(tit!=newTriangles.end())
    {
        visit_points(*tit);
        ++tit;
    }
}

void D25ActiveContours::visit_points( HTriangleSeed &tr )
{
    const float eps=0.001f;

    //Mark the triangle vertices as 'nodes'
    tr.p1->isNode=true;
    tr.p2->isNode=true;
    tr.p3->isNode=true;

    //Select points from the neighborhood. Define the radius
    Vector<float,3> mid=(tr.p1->p+tr.p2->p+tr.p3->p)/3;
    float tmp, rad=maxSurfaceDepth;
    rad=rad>(tmp=float((tr.p1->p-mid).eucl_norm()))?rad:tmp;
    rad=rad>(tmp=float((tr.p2->p-mid).eucl_norm()))?rad:tmp;
    rad=rad>(tmp=float((tr.p3->p-mid).eucl_norm()))?rad:tmp;

    //Pre-search: choose all points in the range
    HPointSeed mids;
    mids.p=mid;
    std::vector<HContainerElement> v;
    vertices->tree.find_within_range(HContainerElement(&mids), rad, std::back_inserter(v));

    //Check if the selected points are inside the prism
    //Calculate triangle prism basis matrix
    Vector<float,3> X=tr.p2->p-tr.p1->p;
    Vector<float,3> Y=tr.p3->p-tr.p1->p;
    Vector<float,3> Z=getNormalVector(tr); Z=Z/float(Z.eucl_norm())*maxSurfaceDepth;
    Vector<float,3> O=tr.p1->p;

    boost::numeric::ublas::matrix<float> m(3,3);
    m(0,0)=X.x(); m(0,1)=Y.x(); m(0,2)=Z.x();
    m(1,0)=X.y(); m(1,1)=Y.y(); m(1,2)=Z.y();
    m(2,0)=X.z(); m(2,1)=Y.z(); m(2,2)=Z.z();
    
    boost::numeric::ublas::matrix<float> im(3,3);  
    
    if( blas::invert_matrix(m,im))
    {
        std::vector<HContainerElement>::iterator it=v.begin();
        while(it!=v.end())
        {
            if(!it->ps->isVisited)
            {
                boost::numeric::ublas::matrix<float> mp(3,1);
                mp(0,0)=it->ps->p.x()-O.x();
                mp(1,0)=it->ps->p.y()-O.y();
                mp(2,0)=it->ps->p.z()-O.z();

                boost::numeric::ublas::matrix<float> am=prod(im,mp);

                float a=am(0,0);
                float b=am(1,0);
                float c=am(2,0);

                if(a+b<=1+eps&&a>=-eps&&b>=-eps&&c>=-1&&c<=1)
                {

                    visit_point(it->ps);	
                }
            }
            ++it;
        }
    }
    
    
}

inline void D25ActiveContours::visit_point( HPointSeed* p )
{
    if(p->isVisited==false)
    {
        p->isVisited=true;
        --unvisitedCount;
    }
}


inline void D25ActiveContours::add_active_edge(const HEdgeSeed &e )
{
    //Init the list of segments
    std::list<HEdgeSeed> segments;
    segments.push_back(e);

    //Check overlapping of the segments with all existing edges
    kill_overlapping_regular_segments(segments,activeEdges);
    kill_overlapping_regular_segments(segments,frozenEdges);

    activeEdges.splice(activeEdges.end(),segments);
}

void D25ActiveContours::kill_overlapping_regular_segments( std::list<HEdgeSeed> &segmentParts, std::list<HEdgeSeed> &edgeList )
{
    std::list<HEdgeSeed> newEdgeList;

    if(segmentParts.size()==0)return;

    std::list<HEdgeSeed>::iterator it=edgeList.begin();
    while(it!=edgeList.end())
    {
        HEdgeSeed e=*it;

        bool isEdgeDeleted=false;

        std::list<HEdgeSeed>::iterator its=segmentParts.begin();
        while(its!=segmentParts.end())
        {
            HEdgeSeed es=*its;

            float t1,t2;
            if(segment_overlap_parameter(*es.p1,e,t1)&&segment_overlap_parameter(*es.p2,e,t2))
            {
                //Set order: t1>=t2
                if(t1<t2)
                {
                    HPointSeed* tmp=(*its).p1;
                    (*its).p1=(*its).p2;
                    (*its).p2=tmp;

                    float ttmp=t1;
                    t1=t2;
                    t2=ttmp;
                }

                if(t1>1)
                {
                    //a--x2--b--x1
                    if(t2>0&&t2<1)
                    {
                        HPointSeed* tmp=(*it).p2;
                        //change edge a--b -> a--x2
                        (*it).p2=(*its).p2;
                        //change segment x2--x1 -> b--x1
                        (*its).p2=tmp;

                        break;
                    }
                    //a,x2--b--x1
                    else if(t2==0)
                    {
                        //change segment x2--x1 -> b--x1
                        (*its).p2=(*it).p2;

                        //delete edge a--b
                        it=edgeList.erase(it);
                        isEdgeDeleted=true;

                        break;
                    }
                    //x2--a--b--x1
                    else if(t2<0)
                    {
                        HEdgeSeed s1,s2;

                        s1.p1=(*its).p1;
                        s1.p2=(*it).p2;
                        s1.propagationVector=(*its).propagationVector;

                        s2.p1=(*it).p1;
                        s2.p2=(*its).p2;
                        s2.propagationVector=(*its).propagationVector;

                        //delete segment x2--x1
                        segmentParts.erase(its);

                        //insert segments x2--a, b--x1
                        segmentParts.push_front(s1);
                        segmentParts.push_front(s2);

                        //delete edge a--b
                        it=edgeList.erase(it);
                        isEdgeDeleted=true;

                        break;
                    }
                }
                else if(t1==1)
                {
                    //a--x2--b,x1
                    if(t2>0&&t2<1)
                    {
                        //change edge a--b -> a--x2
                        (*it).p2=(*its).p2;

                        //delete segment x2--x1
                        segmentParts.erase(its);

                        break;
                    }
                    //a,x2--b,x1
                    else if(t2==0)
                    {
                        //delete segment x2--x1
                        segmentParts.erase(its);

                        //delete edge a--b
                        it=edgeList.erase(it);
                        isEdgeDeleted=true;

                        break;
                    }
                    //x2--a--b,x1
                    else if(t2<0)
                    {
                        //change segment x2--x1 -> x2--a
                        (*its).p1=(*it).p1;

                        //delete edge a--b
                        it=edgeList.erase(it);
                        isEdgeDeleted=true;

                        break;
                    }
                }
                else if(t1>0&&t1<1)
                {
                    //a--x2--x1--b
                    if(t2>0&&t2<1)
                    {
                        HEdgeSeed s1,s2;

                        s1.p1=(*it).p1;
                        s1.p2=(*its).p2;
                        s1.propagationVector=(*it).propagationVector;

                        s2.p1=(*its).p1;
                        s2.p2=(*it).p2;
                        s2.propagationVector=(*it).propagationVector;

                        //delete edge a--b
                        it=edgeList.erase(it);
                        isEdgeDeleted=true;

                        //insert edges a--x2, x1--b
                        newEdgeList.push_back(s1);
                        newEdgeList.push_back(s2);

                        //delete segment x2--x1
                        segmentParts.erase(its);

                        break;
                    }
                    //a,x2--x1--b
                    else if(t2==0)
                    {
                        //change edge a--b -> x1--b
                        (*it).p1=(*its).p1;

                        //delete segment x2--x1
                        segmentParts.erase(its);

                        break;
                    }
                    //x2--a--x1--b
                    else if(t2<0)
                    {
                        HPointSeed* tmp=(*it).p1;
                        //change edge a--b -> x1--b
                        (*it).p1=(*its).p1;
                        //change segment x2--x1 -> x2--a
                        (*its).p1=tmp;

                        break;
                    }
                }
            }

            ++its;
        }

        if(!isEdgeDeleted)++it;
    }

    if(newEdgeList.size()>0)edgeList.splice(edgeList.end(),newEdgeList);

}


bool D25ActiveContours::segment_overlap_parameter(const HPointSeed &ps, const HEdgeSeed &e, float &t)
{
    const float eps=0.001f;

    Vector<float,3> v1=ps.p-e.p1->p;
    Vector<float,3> v2=e.p2->p-e.p1->p;

    if(ps.p==e.p1->p)
    {
        t=0;
        return true;
    }
    if(ps.p==e.p2->p)
    {
        t=1;
        return true;
    }
    if(e.p1->p==e.p2->p)
    {
        return false;
    }

    float dotProduct=v1*v2;
    float normV2=float(v2.eucl_norm());
    float normsProduct=float(v1.eucl_norm())*normV2;

    //If collinear, |cos|~1
    if(dotProduct>=normsProduct-eps||dotProduct<=-normsProduct+eps)
    {
        t=dotProduct/normV2/normV2;
        return true;
    }

    return false;
}

inline bool D25ActiveContours::exclude_small_angles( const HEdgeSeed &e, HPointSeed* &ps )
{

    if(stick_to_adjacent_edge(e, ps, activeEdges) || stick_to_adjacent_edge(e, ps, frozenEdges) )
    {
        return true;
    }

    return false;
}

bool D25ActiveContours::stick_to_adjacent_edge( const HEdgeSeed &e, HPointSeed* &ps, std::list<HEdgeSeed> &edgeList )
{

    std::list<HEdgeSeed>::const_iterator it=edgeList.begin();
    while(it!=edgeList.end())
    {
        HEdgeSeed ee=*it;

        //Excluding the base edge
        if(e.p1->p==ee.p1->p&&e.p2->p==ee.p2->p||
            e.p1->p==ee.p2->p&&e.p2->p==ee.p1->p)
        {
            ++it;
            continue;
        }

        //if at least one of the edges is coinciding, quit
        if((e.p1->p==ee.p1->p||e.p2->p==ee.p1->p)&&ps->p==ee.p2->p||
            (e.p1->p==ee.p2->p||e.p2->p==ee.p2->p)&&ps->p==ee.p1->p)
        {
            return true;
        }

        //Sticking
        Vector<float,3> p;
        if((p=e.p1->p)==ee.p1->p||e.p1->p==ee.p2->p||
            (p=e.p2->p)==ee.p1->p||e.p2->p==ee.p2->p)
        {
            if(p==ee.p2->p)
                ee.swap();
            
            Vector<float,3> v1=ps->p-p;
            double normV1=v1.eucl_norm();

            Vector<float,3> v2=ee.p2->p-p;
            double normV2=v2.eucl_norm();

            if(normV1==0||normV2==0)return false;

            double cosa=v1*v2/float(normV1*normV2);

            if(cosa>maxExcludedAngle)
            {
                ps=ee.p2;
                return true;
            }

        }

        ++it;
    }
    return false;
}


common::Mesh D25ActiveContours::build_mesh(std::vector<Vector<float,3> > &v)
{
    set_vertices(v);

    //Building the set of triangles
    while(grow_step());

    //Generate and return the mesh
    return get_mesh();
}

common::Mesh D25ActiveContours::build_mesh()
{
    //Cleaning all the auxiliary containers
    activeEdges.clear();
    frozenEdges.clear();
    triangles.clear();
    unvisitedCount = static_cast<unsigned>(vertices->tree.size());

    //Building the set of triangles
    while(grow_step());

    //Generate and return the mesh
    return get_mesh();
}



inline bool D25ActiveContours::triangle_mesh_3d_intersection( const HTriangleSeed &t )
{
    Triangle<Vector<float,3> > t1(t.p1->p,t.p2->p,t.p3->p);

    std::list<HTriangleSeed>::iterator tit = triangles.begin();
    while(tit!=triangles.end())
    {
        HTriangleSeed tt=*tit;
        
        Triangle<Vector<float,3> > t2(tt.p1->p,tt.p2->p,tt.p3->p);
        
        if(triangles_3d_intersection(t1,t2))return true;

        ++tit;
    }

    return false;
}

void D25ActiveContours::edge_stitch(HEdgeSeed e )
{
    HPointSeed* pps1=get_propagated_vertex(e,true);

    bool isStitched=false;

    if(pps1)
    {
        Triangle<Vector<float,3> > t1(e.p1->p,e.p2->p,pps1->p);

        for(std::list<HEdgeSeed>::iterator ite=frozenEdges.begin(); ite!=frozenEdges.end(); ++ite)
        {
            HEdgeSeed ee=*ite;

            if(e==ee)continue;

            bool b11=false,b12=false,b21=false,b22=false;

            //If the edge is adjacent
            if((b11=(e.p1==ee.p1))||(b12=(e.p1==ee.p2))||(b21=(e.p2==ee.p1))||(b22=(e.p2==ee.p2)))
            {
                HPointSeed* pps2=get_propagated_vertex(ee,true);
                
                if(pps2)
                {
                    Triangle<Vector<float,3> > t2(ee.p1->p,ee.p2->p,pps2->p);

                    if(triangles_3d_intersection(t1,t2))
                    {
                        //changing all to b11 condition
                        if(b12)ee.swap();
                        else if(b21)e.swap();
                        else if(b22)
                        {
                            e.swap();
                            ee.swap();
                        }

                        HPointSeed *p1,*p2,*p3;
                        p1=e.p1;
                        p2=e.p2;
                        p3=ee.p2;

                        //stitching
                        HEdgeSeed ne(p2,p3);
                        HEdgeSeed nex(p3,p1); //fictive
                        HEdgeSeed nexx(p1,p2); //fictive


                        if(get_edges_propagations(ne,nex,nexx))
                        {
                            HTriangleSeed tr;
                            tr.p1=p1;
                            tr.p2=p2;
                            tr.p3=p3;
 
                            //Test for triangles collisions
                            if(!triangle_mesh_3d_intersection(tr))
                            {
                                //Delete the current frozen neighboring edge
                                frozenEdges.erase(ite);

                                //Add a new active edge
                                add_active_edge(ne);

                                //Add new triangle to the mesh
                                triangles.push_back(tr);

                                isStitched=true;

                                break;
                            }


                        }
                    }

                }

            }
        }
    }


    //Add to frozen edges, if wasn't stitched on this step
    if(!isStitched)
        frozenEdges.push_back(e);

}

Vector<float,3> D25ActiveContours::get_surface_normal( Vector<float,3> p, float windowRadius )
{
    //Select points from the neighborhood
    HPointSeed ps;
    ps.p=p;
    std::vector<HContainerElement> neighbours;
    float const range = windowRadius;
    vertices->tree.find_within_range(HContainerElement(&ps), range, std::back_inserter(neighbours));

    unsigned pointCount = static_cast<unsigned>(neighbours.size());

    //Principal Component Analysis (PCA)
    Vector<float,3> mean(0,0,0);	
    mean=mean/(float)pointCount;

    double** A=new double*[6]; 
    for(int i=0; i<6; ++i)
    {
        A[i]=new double[3];
        A[i][0]=A[i][1]=A[i][2]=0;
    }

    double* S2=new double[3];

    std::vector<HContainerElement>::const_iterator itp=neighbours.begin();
    while(itp!=neighbours.end())
    {
        Vector<float,3> pp=(*itp).ps->p;

        A[0][0]+=(pp.x()-mean.x())*(pp.x()-mean.x()); A[1][0]+=(pp.y()-mean.y())*(pp.x()-mean.x()); A[2][0]+=(pp.z()-mean.z())*(pp.x()-mean.x()); 
        A[0][1]+=(pp.x()-mean.x())*(pp.y()-mean.y()); A[1][1]+=(pp.y()-mean.y())*(pp.y()-mean.y()); A[2][1]+=(pp.z()-mean.z())*(pp.y()-mean.y());
        A[0][2]+=(pp.x()-mean.x())*(pp.z()-mean.z()); A[1][2]+=(pp.y()-mean.y())*(pp.z()-mean.z()); A[2][2]+=(pp.z()-mean.z())*(pp.z()-mean.z());

        ++itp;
    }
    for(int i=0; i<3; ++i)
    {
        A[i][0]/=pointCount;
        A[i][1]/=pointCount;
        A[i][2]/=pointCount;
    }

    int n=3;

    //Calculate eigenvectors using SVD decomposition
    svd(A,S2,n);

    int maxIndex=(S2[0]<S2[1])?(S2[0]<S2[2]?0:2):(S2[1]<S2[2]?1:2);

    Vector<float,3> v((float)A[3+maxIndex][0], (float)A[3+maxIndex][1], (float)A[3+maxIndex][2]);

    for(int i=0; i<6; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] S2;

    return v/float(v.eucl_norm());
}



const std::vector<HPointSeed>* D25ActiveContours::get_vertices()
{
    return &vertices->linear;
}

const std::list<HEdgeSeed>* D25ActiveContours::get_active_edges()
{
    return &activeEdges;
}

const std::list<HEdgeSeed>* D25ActiveContours::get_frozen_edges()
{
    return &frozenEdges;
}

const std::list<HTriangleSeed>* D25ActiveContours::get_triangles()
{
    return &triangles;
}



bool D25ActiveContours::triangles_3d_intersection(const Triangle<Vector<float,3> > &t1,
                                                  const Triangle<Vector<float,3> > &t2)
{
    const float eps=0.001f;
    float alpha=tetrahedronBaseAngle<eps?eps:tetrahedronBaseAngle;

    //if the triangles are adjacent
    if(t1.A()==t2.A() || t1.A()==t2.B() || t1.A()==t2.C() ||
        t1.B()==t2.A() || t1.B()==t2.B() || t1.B()==t2.C() ||
        t1.C()==t2.A() || t1.C()==t2.B() || t1.C()==t2.C())
    {

        TriangularDipyramid tdp1=TriangularDipyramid::from_triangle_and_angle(t1,alpha);
        TriangularDipyramid tdp2=TriangularDipyramid::from_triangle_and_angle(t2,alpha);

        return tdp1.intersects(tdp2);
    }
    else
    {
        TriangularDipyramid tdp1=TriangularDipyramid::from_triangle_and_height(t1,maxSurfaceDepth);
        TriangularDipyramid tdp2=TriangularDipyramid::from_triangle_and_height(t2,maxSurfaceDepth);

        return tdp1.intersects(tdp2);
    }
}

bool D25ActiveContours::triangle_degenerate( const HTriangleSeed &t )
{
    Vector<float,3> v1a=t.p2->p-t.p1->p;
    Vector<float,3> v1b=t.p3->p-t.p1->p;

    Vector<float,3> v2a=t.p1->p-t.p2->p;
    Vector<float,3> v2b=t.p3->p-t.p2->p; 

    Vector<float,3> v3a=t.p1->p-t.p3->p;
    Vector<float,3> v3b=t.p2->p-t.p3->p;
 
    float cos1=v1a*v1b/float((v1a.eucl_norm()*v1b.eucl_norm()));
    float cos2=v2a*v2b/float((v2a.eucl_norm()*v2b.eucl_norm()));
    float cos3=v3a*v3b/float((v3a.eucl_norm()*v3b.eucl_norm()));

    if(fabs(cos1)>maxExcludedAngle||fabs(cos2)>maxExcludedAngle||fabs(cos3)>maxExcludedAngle)return true;
    else return false;
}

void D25ActiveContours::set_vertices( std::vector<Vector<float,3> > &v )
{
    //Cleaning all the auxiliary containers
    activeEdges.clear();
    frozenEdges.clear();
    triangles.clear();

    //Create and fill in a vertex container
    if(vertices)delete vertices;
    vertices=new HVertexContainer(v);

    unvisitedCount = static_cast<unsigned>(vertices->tree.size());
}

bool D25ActiveContours::grow_step()
{	
    if(unvisitedCount==0&&activeEdges.size()==0)
    {
        unsigned int frozenBefore = static_cast<unsigned>(frozenEdges.size());
        post_stitch();
        return activeEdges.size()>0 || frozenEdges.size()<frozenBefore;
    }
    else if(activeEdges.size()>0)model_grow();
    else model_init();

    return true;
}

void D25ActiveContours::post_stitch()
{
    if(frozenEdges.size()==0)return;

    for(std::list<HEdgeSeed>::iterator it=frozenEdges.begin(); it!=frozenEdges.end(); ++it )
    {
        HEdgeSeed e=*it;
        HPointSeed* pps1=get_propagated_vertex(e,true);
        if(!pps1)continue;

        bool isStitched=false;

        for(std::list<HEdgeSeed>::iterator ite=frozenEdges.begin(); ite!=frozenEdges.end(); ++ite )
        {
            HEdgeSeed ee=*ite;

            if(e==ee)continue;
            HPointSeed* pps2=get_propagated_vertex(ee,true);
            if(!pps2)continue;

            bool b11=false,b12=false,b21=false,b22=false;

            //If the edge is adjacent
            if((b11=(e.p1==ee.p1))||(b12=(e.p1==ee.p2))||(b21=(e.p2==ee.p1))||(b22=(e.p2==ee.p2)))
            {
                //changing all to b11 condition
                if(b12)ee.swap();
                else if(b21)e.swap();
                else if(b22)
                {
                    e.swap();
                    ee.swap();
                }

                //Check the stitching condition (see the "Red Notebook")
                {
                    Vector<float,3> a=e.p2->p-e.p1->p;
                    Vector<float,3> b=ee.p2->p-ee.p1->p;
                    Vector<float,3> norm=getNormalVector(a,b);

                    Vector<float,3> n1=getNormalVector(a,norm);
                    float signCor1=n1*e.propagationVector>0?1.0f:-1.0f;
                    n1=n1*signCor1;

                    Vector<float,3> n2=getNormalVector(b,norm);
                    float signCor2=n2*ee.propagationVector>0?1.0f:-1.0f;
                    n2=n2*signCor2;

                    float sgn1=a*n2;
                    float sgn2=b*n1;
                    float cosab=a*b/float((a.eucl_norm()*b.eucl_norm()));

                    if(sgn1>0&&sgn2>0)
                        if(cosab>maxStitchedAngle)
                        {
                            HTriangleSeed tr;
                            tr.p1=e.p1;
                            tr.p2=e.p2;
                            tr.p3=ee.p2;

                            HEdgeSeed ne;
                            ne.p1=e.p2;
                            ne.p2=ee.p2;

                            //Delete the current frozen neighboring edge
                            frozenEdges.erase(ite);

                            //Add a new active edge
                            add_active_edge(ne);

                            //Add new triangle to the mesh
                            triangles.push_back(tr);

                            isStitched=true;
                            break;

                        }

                }

            }
        }

        if(isStitched)
        {
            frozenEdges.remove(e);
            break;
        }

    }

} // namespace surfaces

common::Mesh D25ActiveContours::get_mesh()
{
    //Construct a mesh
    common::Mesh m(triangles.size());

    //Create a reference map (from the local triangles nodes to the mesh vertices)
    std::map<HPointSeed*,size_t> mymap;
    std::list<HTriangleSeed>::const_iterator itt=triangles.begin();
    while(itt!=triangles.end())
    {
        for(int j=0; j<3; ++j)
        {
            HPointSeed* ps=(j==0)?(itt->p1):(j==1?itt->p2:itt->p3);

            //If the triangle node is not yet in the map
            if(mymap.find(ps)==mymap.end())
            {	
                //Add the vertex/node into the mesh and to the reference map
                size_t ind=m.add_vertex(common::Mesh::Vertex(ps->p.x(),ps->p.y(),ps->p.z()));
                mymap[ps]=ind;
            }
        }

        //Add the face from the processed triangle nodes
        size_t A=mymap[itt->p1];
        size_t B=mymap[itt->p2];
        size_t C=mymap[itt->p3];

        m.add_face(common::Mesh::Face(A,B,C));

        ++itt;
    }

    return m;
}


} // namespace surfaces

} // namespace methods
} // namespace common
