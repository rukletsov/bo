#include "stdafx.h"

#include "d25_active_contours.hpp"

#include <functional>

#include "svd.h"
#include "kdtree++/kdtree.hpp"


//Uncomment for equal-square triangles propagation  
//#define USE_EQUAL_SQUARE_PROPAGATION

/*! \struct HTreeElement
\brief HPointSeed wrapper for HVertexContainer utilization 
*/
struct HContainerElement
{
	/*! Default constructor */
	HContainerElement():ps(0){}

	/*! Constructor */
	HContainerElement(HPointSeed* ps):ps(ps){}

	/*! Access operator */
	inline float operator [] (const int t) const 
	{
		return (*ps)[t];
	}

	/*! Comparison operator */
	inline bool operator == (const HContainerElement &other) const
	{
		return (*ps)==(*other.ps);
	}

	/*! Pointer to an elementary vertex item*/
	HPointSeed* ps;
};


//Brackets accessor
inline float bac( HContainerElement t, size_t k ) { return t[k]; }

//3D Tree type
typedef KDTree::KDTree<3, HContainerElement, std::pointer_to_binary_function<HContainerElement,size_t,float>> D3Tree;

//3D Tree Wrapper
class HVertexContainer
{
public:

	HVertexContainer(std::list<HPoint3<float>>& vertices):tree(D3Tree(std::ptr_fun(bac)))
	{		
		linear.resize(vertices.size());

		//Filling in the 3D Tree
		int cnt=0;
		for(std::list<HPoint3<float>>::const_iterator itp=vertices.begin(); itp!=vertices.end(); ++itp)
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
	std::vector<HPointSeed> linear;
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
		return ((!ce.ps->isVisited)||(checkNodes&&ce.ps->isNode)||(checkVisited&&ce.ps->isVisited))&&((searchCenter.p-ce.ps->p).getNorm()>minDistance);
	}
protected:
	HPointSeed searchCenter;
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
		ba=ce2.ps->p-ce1.ps->p;
		absBa=ba.getNorm();
	}
	inline bool operator()( HContainerElement const& ce ) const
	{
		bool isPretender=(!ce.ps->isVisited)||(checkNodes&&ce.ps->isNode)||(checkVisited&&ce.ps->isVisited);
		if(!isPretender)return false;

		if(absBa==0)return false;

		//Collinearity test
		HPoint3<float> ca=ce.ps->p-ce1.ps->p;;
		float absCa2=ca.x*ca.x+ca.y*ca.y+ca.z*ca.z;

		float scalBaCa=ba*ca;
		float projCaOnBa=scalBaCa/absBa;

		float residual2=absCa2-projCaOnBa*projCaOnBa;

		if(residual2>eps)return true;
		else return false;
	}
protected:
	HContainerElement ce1;
	HContainerElement ce2;
	bool checkNodes;
	bool checkVisited;
	float eps;
	HPoint3<float> ba;
	float absBa;
};


HActiveContours::HActiveContours()
{
	//Default init

	minInitDistance=8.0f;
	maxInitDistance=1.5f*minInitDistance;
	maxProjectionNodeDistance=0.5f*minInitDistance;
	normalNeighborhoodRadius=maxInitDistance;
	maxSurfaceDepth=5.0f;

	maxExcludedAngle=0.85f;
	maxStitchedAngle=-0.90f;

	faceSurfaceFactor=0.5;

	vertices=0;
}

HActiveContours::~HActiveContours()
{
	delete vertices;
}

inline HPointSeed* HActiveContours::getClosestPoint( const HPointSeed &ps, bool checkNodes, bool checkVisited )
{
	HContainerElement ce(0);

	PredicateClosestPointWithMinDistance pred(HContainerElement(const_cast<HPointSeed*>(&ps)),minInitDistance,checkNodes,checkVisited);
	std::pair<D3Tree::const_iterator,float> nif = vertices->tree.find_nearest_if(HContainerElement(const_cast<HPointSeed*>(&ps)),maxInitDistance,pred);
	if(nif.first!=vertices->tree.end())ce=*nif.first;

	return ce.ps;
}

HPointSeed* HActiveContours::getClosestMinFuncPoint( const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited )
{
	HPoint3<float> v1=ps2.p-ps1.p;
	float a=v1.getNorm();

	HPointSeed mid;
	mid.p=(ps1.p+ps2.p)/2;

	//Pre-search: choose all points in the range
	std::vector<HContainerElement> v;
	float const range = maxInitDistance;
	vertices->tree.find_within_range(HContainerElement(&mid), range, std::back_inserter(v));

	HContainerElement ce(0);
	float func=-1;

	std::vector<HContainerElement>::const_iterator it = v.begin();
	while(it!=v.end())
	{
		if((!(*it).ps->isVisited)||(checkNodes&&(*it).ps->isNode)||(checkVisited&&(*it).ps->isVisited))
		{	
			HPoint3<float> v2=(*it).ps->p-ps2.p;
			HPoint3<float> v3=ps1.p-(*it).ps->p;

			float b=v2.getNorm();
			float c=v3.getNorm();

			if(b<maxInitDistance&&c<maxInitDistance&&b>minInitDistance&&c>minInitDistance)
			{
				float f=abs(a-b)+abs(a-c);

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


inline HPointSeed* HActiveContours::getClosestNoncollinearPoint(const HPointSeed &ps, const HPointSeed &ps1, const HPointSeed& ps2, bool checkNodes, bool checkVisited)
{
	HContainerElement ce(0);

	PredicateClosestPointNonCollinear pred(HContainerElement(const_cast<HPointSeed*>(&ps1)),HContainerElement(const_cast<HPointSeed*>(&ps2)),checkNodes,checkVisited);
	std::pair<D3Tree::const_iterator,float> nif = vertices->tree.find_nearest_if(HContainerElement(const_cast<HPointSeed*>(&ps)),maxProjectionNodeDistance,pred);
	if(nif.first!=vertices->tree.end())ce=*nif.first;

	return ce.ps;
}


inline float HActiveContours::getDistance( const HPointSeed &ps1, const HPointSeed &ps2 )
{
	return (ps1.p-ps2.p).getNorm();
}


void HActiveContours::modelInit()
{
	HPointSeed *pps1=0,*pps2=0,*pps3=0;

	D3Tree::mutable_iterator it = vertices->tree.begin();
	while(it!=vertices->tree.end())
	{
		if(!(*it).ps->isVisited)
		{
			pps1=(*it).ps;

			visitPoint(pps1);

			pps2=getClosestPoint(*pps1, true, false);

			if(pps2)
			{
				visitPoint(pps2);

				pps3=getClosestMinFuncPoint(*pps1, *pps2, true, false);

				if(pps3)
				{
					HTriangleSeed tr;
					tr.p1=pps1;
					tr.p2=pps2;
					tr.p3=pps3;

					visitPoint(pps3);

					if(!triangleDegenerate(tr)&&!triangleMesh3DIntersection(tr))
					{
						HEdgeSeed e1, e2, e3;

						e1.p1=e3.p2=pps1;
						e2.p1=e1.p2=pps2;
						e3.p1=e2.p2=pps3;

						if(getEdgesPropagations(e1,e2,e3))
						{
							addActiveEdge(e1);
							addActiveEdge(e2);
							addActiveEdge(e3);

							visitPoints(tr);

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


void HActiveContours::modelGrow()
{
	if(activeEdges.size()==0)return;

	HEdgeSeed e=*activeEdges.begin();

	HPointSeed* pps=getPropagatedVertex(e, false);

	if(pps)
	{
		//Exclude small angles
		excludeSmallAngles(e, pps);

		HTriangleSeed tr;
		tr.p1=e.p1;
		tr.p2=e.p2;
		tr.p3=pps;

		if(triangleDegenerate(tr))
		{
			frozenEdges.push_back(e);
		}
		//Test for triangles collisions
		else if(!triangleMesh3DIntersection(tr))
		{
			//New edges based on an old one and the propagated point
			HEdgeSeed e1(e.p2,pps);
			HEdgeSeed e2(pps,e.p1);
			HEdgeSeed e3=e; //fictive

			//If the propagations are successfully calculated
			if(getEdgesPropagations(e1,e2,e3))
			{
				//Add two new active edges
				addActiveEdge(e1);
				addActiveEdge(e2);

				//Mark visited and under projection points
				visitPoints(tr);	

				//Add new triangle to the mesh
				triangles.push_back(tr);
			}

		}
		else 
		{
			//"On-the-fly" stitching
			edgeStitch(e);
		}
	}
	else 
	{
		//"On-the-fly" stitching
		edgeStitch(e);		
	}

	//Delete the processed edge from the active edges if it is yet here
	std::list<HEdgeSeed>::iterator it=activeEdges.begin();
	if(it!=activeEdges.end()&&e==*it)
		activeEdges.erase(it);
}

bool HActiveContours::getEdgesPropagations( HEdgeSeed &e1, HEdgeSeed &e2, HEdgeSeed &e3 )
{
	if(e1.p1!=e3.p2 || e2.p1!=e1.p2 || e3.p1!=e2.p2)return false;

	HPoint3<float> mid=(e1.p1->p+e2.p1->p+e3.p1->p)/3;

	HPoint3<float> v1=e1.p2->p-e1.p1->p;
	HPoint3<float> v2=e2.p2->p-e2.p1->p;
	HPoint3<float> v3=e3.p2->p-e3.p1->p;

	HPoint3<float> propagationE1,propagationE2,propagationE3;

	//PCA-based approximation 
	if(faceSurfaceFactor!=0)
	{
		//PCA-based surface normal calculation
		HPoint3<float> midNorm=getSurfaceNormal(mid,normalNeighborhoodRadius);

		//Propagation directions. Surface normal calculation
		propagationE1=getNormalVector(v1,midNorm);
		propagationE2=getNormalVector(v2,midNorm);
		propagationE3=getNormalVector(v3,midNorm);

		//Outer directions correction
		HPoint3<float> medianE1=(e1.p1->p+e1.p2->p)/2-mid;
		HPoint3<float> medianE2=(e2.p1->p+e2.p2->p)/2-mid;
		HPoint3<float> medianE3=(e3.p1->p+e3.p2->p)/2-mid;

		float cosMP1=(medianE1*propagationE1)/(medianE1.getNorm()*propagationE1.getNorm());
		float cosMP2=(medianE2*propagationE2)/(medianE2.getNorm()*propagationE2.getNorm());
		float cosMP3=(medianE3*propagationE3)/(medianE3.getNorm()*propagationE3.getNorm());

		if(cosMP1<0)propagationE1=propagationE1*(-1);
		if(cosMP2<0)propagationE2=propagationE2*(-1);
		if(cosMP3<0)propagationE3=propagationE3*(-1);

		//Mixture with the PCA-face propagations
		propagationE1=(propagationE1/propagationE1.getNorm())*faceSurfaceFactor+(medianE1/medianE1.getNorm())*(1-faceSurfaceFactor);
		propagationE2=(propagationE2/propagationE2.getNorm())*faceSurfaceFactor+(medianE2/medianE2.getNorm())*(1-faceSurfaceFactor);
		propagationE3=(propagationE3/propagationE3.getNorm())*faceSurfaceFactor+(medianE3/medianE3.getNorm())*(1-faceSurfaceFactor);
	}
	//Triangle face based approximation
	else
	{
		propagationE1=(e1.p1->p+e1.p2->p)/2-mid;
		propagationE2=(e2.p1->p+e2.p2->p)/2-mid;
		propagationE3=(e3.p1->p+e3.p2->p)/2-mid;
	}


	//Normalization
	propagationE1=propagationE1/propagationE1.getNorm();
	propagationE2=propagationE2/propagationE2.getNorm();
	propagationE3=propagationE3/propagationE3.getNorm();

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


inline HPointSeed* HActiveContours::getPropagatedVertex(const HEdgeSeed &e, bool checkVisited)
{
	HPointSeed p;

	HPoint3<float> mid=(e.p1->p+e.p2->p)/2;

	p.p=mid+e.propagationVector;

	HPointSeed* ps=getClosestNoncollinearPoint(p,*e.p1,*e.p2,true, checkVisited);

	return ps;
}


void HActiveContours::visitPoints( std::list<HTriangleSeed> &newTriangles )
{
	std::list<HTriangleSeed>::iterator tit = newTriangles.begin();
	while(tit!=newTriangles.end())
	{
		visitPoints(*tit);
		++tit;
	}
}

void HActiveContours::visitPoints( HTriangleSeed &tr )
{
	const float eps=0.001f;

	//Mark the triangle vertices as 'nodes'
	tr.p1->isNode=true;
	tr.p2->isNode=true;
	tr.p3->isNode=true;

	//Select points from the neighborhood. Define the radius
	HPoint3<float> mid=(tr.p1->p+tr.p2->p+tr.p3->p)/3;
	float tmp, rad=maxSurfaceDepth;
	rad=rad>(tmp=(tr.p1->p-mid).getNorm())?rad:tmp;
	rad=rad>(tmp=(tr.p2->p-mid).getNorm())?rad:tmp;
	rad=rad>(tmp=(tr.p3->p-mid).getNorm())?rad:tmp;

	//Pre-search: choose all points in the range
	HPointSeed mids;
	mids.p=mid;
	std::vector<HContainerElement> v;
	vertices->tree.find_within_range(HContainerElement(&mids), rad, std::back_inserter(v));

	//Check if the selected points are inside the prism
	//Calculate triangle prism basis matrix
	HPoint3<float> X=tr.p2->p-tr.p1->p;
	HPoint3<float> Y=tr.p3->p-tr.p1->p;
	HPoint3<float> Z=getNormalVector(tr); Z=Z/Z.getNorm()*maxSurfaceDepth;
	HPoint3<float> O=tr.p1->p;

	MatrixTCL<float> m(3,3);
	m(0,0)=X.x; m(0,1)=Y.x; m(0,2)=Z.x;
	m(1,0)=X.y; m(1,1)=Y.y; m(1,2)=Z.y;
	m(2,0)=X.z; m(2,1)=Y.z; m(2,2)=Z.z;
	m=m.Inv();

	std::vector<HContainerElement>::iterator it=v.begin();
	while(it!=v.end())
	{
		if(!it->ps->isVisited)
		{
			MatrixTCL<float> mp(3,1);
			mp(0,0)=it->ps->p.x-O.x;
			mp(1,0)=it->ps->p.y-O.y;
			mp(2,0)=it->ps->p.z-O.z;

			MatrixTCL<float> am=m*mp;

			float a=am(0,0);
			float b=am(1,0);
			float c=am(2,0);

			if(a+b<=1+eps&&a>=-eps&&b>=-eps&&c>=-1&&c<=1)
			{

				visitPoint(it->ps);	
			}
		}
		++it;
	}
}

inline void HActiveContours::visitPoint( HPointSeed* p )
{
	if(p->isVisited==false)
	{
		p->isVisited=true;
		--unvisitedCount;
	}
}


inline HPoint3<float> HActiveContours::getNormalVector(const HTriangleSeed &t)
{
	HPoint3<float> A=t.p2->p-t.p1->p;
	HPoint3<float> B=t.p3->p-t.p1->p;

	HPoint3<float> p;
	p.x=A.y*B.z-B.y*A.z;
	p.y=A.z*B.x-B.z*A.x;
	p.z=A.x*B.y-B.x*A.y;

	return p;
}

inline void HActiveContours::addActiveEdge(const HEdgeSeed &e )
{
	//Init the list of segments
	std::list<HEdgeSeed> segments;
	segments.push_back(e);

	//Check overlapping of the segments with all existing edges
	killOverlappingRegularSegments(segments,activeEdges);
	killOverlappingRegularSegments(segments,frozenEdges);

	activeEdges.splice(activeEdges.end(),segments);
}

void HActiveContours::killOverlappingRegularSegments( std::list<HEdgeSeed> &segmentParts, std::list<HEdgeSeed> &edgeList )
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
			if(segmentOverlapParameter(*es.p1,e,t1)&&segmentOverlapParameter(*es.p2,e,t2))
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


bool HActiveContours::segmentOverlapParameter(const HPointSeed &ps, const HEdgeSeed &e, float &t)
{
	const float eps=0.001f;

	HPoint3<float> v1=ps.p-e.p1->p;
	HPoint3<float> v2=e.p2->p-e.p1->p;

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
	float normV2=v2.getNorm();
	float normsProduct=v1.getNorm()*normV2;

	//If collinear, |cos|~1
	if(dotProduct>=normsProduct-eps||dotProduct<=-normsProduct+eps)
	{
		t=dotProduct/normV2/normV2;
		return true;
	}

	return false;
}

inline bool HActiveContours::excludeSmallAngles( const HEdgeSeed &e, HPointSeed* &ps )
{

	if(stickToAdjacentEdge(e, ps, activeEdges) || stickToAdjacentEdge(e, ps, frozenEdges) )
	{
		return true;
	}

	return false;
}

bool HActiveContours::stickToAdjacentEdge( const HEdgeSeed &e, HPointSeed* &ps, std::list<HEdgeSeed> &edgeList )
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
		HPoint3<float> p;
		if((p=e.p1->p)==ee.p1->p||e.p1->p==ee.p2->p||
			(p=e.p2->p)==ee.p1->p||e.p2->p==ee.p2->p)
		{
			if(p==ee.p2->p)
			{
				HPointSeed *tmp=ee.p1;
				ee.p1=ee.p2;
				ee.p2=tmp;
			}

			HPoint3<float> v1=ps->p-p;
			float normV1=v1.getNorm();

			HPoint3<float> v2=ee.p2->p-p;
			float normV2=v2.getNorm();

			if(normV1==0||normV2==0)return false;

			float cosa=(v1*v2/normV1/normV2);

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


std::list<HTriangle3<float>> HActiveContours::buildMesh(std::list<HPoint3<float>> &vertexList)
{
	if(vertices)delete vertices;
	vertices=new HVertexContainer(vertexList);

	//Cleaning all the auxiliary containers
	activeEdges.clear();
	frozenEdges.clear();
	triangles.clear();

	unvisitedCount=vertices->tree.size();

	//Building the mesh
	while(growStep());


	//Return the mesh
	std::list<HTriangle3<float>> mesh;
	std::list<HTriangleSeed>::const_iterator itt=triangles.begin();
	while(itt!=triangles.end())
	{
		HTriangle3<float> tr;
		tr.A=itt->p1->p;
		tr.B=itt->p2->p;
		tr.C=itt->p3->p;
		mesh.push_back(tr);
		++itt;
	}
	return mesh;

}

inline bool HActiveContours::triangleMesh3DIntersection( const HTriangleSeed &t )
{
	HTriangle3<float> t1,t2;

	t1.A=t.p1->p;
	t1.B=t.p2->p;
	t1.C=t.p3->p;

	std::list<HTriangleSeed>::iterator tit = triangles.begin();
	while(tit!=triangles.end())
	{
		HTriangleSeed tt=*tit;

		t2.A=tt.p1->p;
		t2.B=tt.p2->p;
		t2.C=tt.p3->p;

		if(triangles3DIntersection(t1,t2))return true;

		++tit;
	}

	return false;
}


inline HPoint3<float> HActiveContours::getNormalVector( const HTriangle3<float> &t )
{
	HPoint3<float> A=t.B-t.A;
	HPoint3<float> B=t.C-t.A;

	return getNormalVector(A,B);
}

inline HPoint3<float> HActiveContours::getNormalVector( const HPoint3<float> A, const HPoint3<float> B )
{
	HPoint3<float> p;

	p.x=A.y*B.z-B.y*A.z;
	p.y=A.z*B.x-B.z*A.x;
	p.z=A.x*B.y-B.x*A.y;

	return p;
}

void HActiveContours::edgeStitch(HEdgeSeed e )
{
	HPointSeed* pps1=getPropagatedVertex(e,true);

	bool isStitched=false;

	if(pps1)
	{
		HTriangle3<float> t1;

		t1.A.x=e.p1->p.x; t1.B.x=e.p2->p.x; t1.C.x=pps1->p.x;
		t1.A.y=e.p1->p.y; t1.B.y=e.p2->p.y; t1.C.y=pps1->p.y;
		t1.A.z=e.p1->p.z; t1.B.z=e.p2->p.z; t1.C.z=pps1->p.z;

		for(std::list<HEdgeSeed>::iterator ite=frozenEdges.begin(); ite!=frozenEdges.end(); ++ite)
		{
			HEdgeSeed ee=*ite;

			if(e==ee)continue;

			bool b11=false,b12=false,b21=false,b22=false;

			//If the edge is adjacent
			if((b11=(e.p1==ee.p1))||(b12=(e.p1==ee.p2))||(b21=(e.p2==ee.p1))||(b22=(e.p2==ee.p2)))
			{
				HPointSeed* pps2=getPropagatedVertex(ee,true);

				if(pps2)
				{
					HTriangle3<float> t2;

					t2.A.x=ee.p1->p.x; t2.B.x=ee.p2->p.x; t2.C.x=pps2->p.x;
					t2.A.y=ee.p1->p.y; t2.B.y=ee.p2->p.y; t2.C.y=pps2->p.y;
					t2.A.z=ee.p1->p.z; t2.B.z=ee.p2->p.z; t2.C.z=pps2->p.z;

					if(triangles3DIntersection(t1,t2))
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


						if(getEdgesPropagations(ne,nex,nexx))
						{
							HTriangleSeed tr;
							tr.p1=p1;
							tr.p2=p2;
							tr.p3=p3;

							//Test for triangles collisions
							if(!triangleMesh3DIntersection(tr))
							{
								//Delete the current frozen neighboring edge
								frozenEdges.erase(ite);

								//Add a new active edge
								addActiveEdge(ne);

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

HPoint3<float> HActiveContours::getSurfaceNormal( HPoint3<float> p, float windowRadius )
{
	//Select points from the neighborhood
	HPointSeed ps;
	ps.p=p;
	std::vector<HContainerElement> neighbours;
	float const range = windowRadius;
	vertices->tree.find_within_range(HContainerElement(&ps), range, std::back_inserter(neighbours));

	int pointCount=neighbours.size();

	//Principal Component Analysis (PCA)
	HPoint3<float> mean(0,0,0);	
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
		HPoint3<float> pp=(*itp).ps->p;

		A[0][0]+=(pp.x-mean.x)*(pp.x-mean.x); A[1][0]+=(pp.y-mean.y)*(pp.x-mean.x); A[2][0]+=(pp.z-mean.z)*(pp.x-mean.x); 
		A[0][1]+=(pp.x-mean.x)*(pp.y-mean.y); A[1][1]+=(pp.y-mean.y)*(pp.y-mean.y); A[2][1]+=(pp.z-mean.z)*(pp.y-mean.y);
		A[0][2]+=(pp.x-mean.x)*(pp.z-mean.z); A[1][2]+=(pp.y-mean.y)*(pp.z-mean.z); A[2][2]+=(pp.z-mean.z)*(pp.z-mean.z);

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

	HPoint3<float> v((float)A[3+maxIndex][0], (float)A[3+maxIndex][1], (float)A[3+maxIndex][2]);

	for(int i=0; i<6; ++i)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] S2;

	return v/v.getNorm();
}

void HActiveContours::setMinInitDistance( float minInitDistance )
{
	this->minInitDistance=minInitDistance;
}

void HActiveContours::setMaxInitDistance( float maxInitDistance )
{
	this->maxInitDistance=maxInitDistance;
}

void HActiveContours::setMaxProjectionNodeDistance( float maxProjectionNodeDistance )
{
	this->maxProjectionNodeDistance=maxProjectionNodeDistance;
}

void HActiveContours::setMaxSurfaceDepth( float maxSurfaceDepth )
{
	this->maxSurfaceDepth=maxSurfaceDepth;
}

void HActiveContours::setMaxExcludedAngle( float maxExcludedAngle )
{
	this->maxExcludedAngle=maxExcludedAngle;
}

void HActiveContours::setNormalNeighborhoodRadius( float normalNeighborhoodRadius )
{
	this->normalNeighborhoodRadius=normalNeighborhoodRadius;
}

void HActiveContours::setFaceSurfaceFactor( float faceSurfaceFactor )
{
	this->faceSurfaceFactor=faceSurfaceFactor;
}

void HActiveContours::setMaxStitchedAngle( float maxStitchedAngle )
{
	this->maxStitchedAngle=maxStitchedAngle;
}

const std::vector<HPointSeed>* HActiveContours::getVerticesVector()
{
	return &vertices->linear;
}

const std::list<HEdgeSeed>* HActiveContours::getActiveEdgeList()
{
	return &activeEdges;
}

const std::list<HEdgeSeed>* HActiveContours::getFrozenEdgeList()
{
	return &frozenEdges;
}

const std::list<HTriangleSeed>* HActiveContours::getTriangleList()
{
	return &triangles;
}

void HActiveContours::initMeshScale( float minInitDistance )
{
	this->minInitDistance=minInitDistance;
	maxInitDistance=1.5f*minInitDistance;
	maxProjectionNodeDistance=0.5f*minInitDistance;
	normalNeighborhoodRadius=maxInitDistance;
}

bool HActiveContours::triangles3DIntersection( const HTriangle3<float> &t1, const HTriangle3<float> &t2 )
{
	const float eps=0.001f;

	HPoint3<float> vertices[2][5];
	HTriangle3<float> faces[2][6];
	HPoint3<float> mid[2];

	//Pyramids construction
	HPoint3<float> z1=getNormalVector(t1);
	HPoint3<float> z2=getNormalVector(t2);

	z1=z1/z1.getNorm()*maxSurfaceDepth/2;
	z2=z2/z2.getNorm()*maxSurfaceDepth/2;

	mid[0]=(t1.A+t1.B+t1.C)/3.0;
	HPoint3<float> p11=mid[0]+z1;
	HPoint3<float> p12=mid[0]-z1;
	mid[1]=(t2.A+t2.B+t2.C)/3.0;
	HPoint3<float> p21=mid[1]+z2;
	HPoint3<float> p22=mid[1]-z2;

	vertices[0][0]=t1.A; vertices[0][1]=t1.B; vertices[0][2]=t1.C;
	vertices[0][3]=p11; vertices[0][4]=p12;
	faces[0][0]=HTriangle3<float>(t1.A, t1.B, p11);
	faces[0][1]=HTriangle3<float>(t1.B, t1.C, p11);
	faces[0][2]=HTriangle3<float>(t1.C, t1.A, p11);
	faces[0][3]=HTriangle3<float>(t1.A, t1.B, p12);
	faces[0][4]=HTriangle3<float>(t1.B, t1.C, p12);
	faces[0][5]=HTriangle3<float>(t1.C, t1.A, p12);

	vertices[1][0]=t2.A; vertices[1][1]=t2.B; vertices[1][2]=t2.C;
	vertices[1][3]=p21; vertices[1][4]=p22;
	faces[1][0]=HTriangle3<float>(t2.A, t2.B, p21);
	faces[1][1]=HTriangle3<float>(t2.B, t2.C, p21);
	faces[1][2]=HTriangle3<float>(t2.C, t2.A, p21);
	faces[1][3]=HTriangle3<float>(t2.A, t2.B, p22);
	faces[1][4]=HTriangle3<float>(t2.B, t2.C, p22);
	faces[1][5]=HTriangle3<float>(t2.C, t2.A, p22);

	for(unsigned i=0; i<2; ++i)
		for(unsigned int t=0; t<6; ++t)
		{
			HTriangle3<float> face=faces[i][t];
			HPoint3<float> norm=getNormalVector(face);
			HPoint3<float> homeVector=mid[i]-face.A;
			float homeDot=(norm*homeVector);
			int homeSign=homeDot>eps?1:(homeDot<-eps?-1:0);

			bool isSeparation=true;
			for(int k=0; k<5; ++k)
			{
				HPoint3<float> alienVector=vertices[1-i][k]-face.A;
				float alienDot=norm*alienVector;
				int alienSign=alienDot>eps?1:(alienDot<-eps?-1:0);
				if(alienSign==homeSign&&alienSign!=0)
				{
					isSeparation=false;
					break;
				}
			}

			if(isSeparation)return false;
		}

		return true;
}

bool HActiveContours::triangleDegenerate( const HTriangleSeed &t )
{
	const float eps=0.001f;
	HPoint3<float> v1=t.p2->p-t.p1->p;
	HPoint3<float> v2=t.p3->p-t.p2->p;
	HPoint3<float> v3=t.p1->p-t.p3->p;

	float cos1=v1*v2/(v1.getNorm()*v2.getNorm());
	float cos2=v2*v3/(v2.getNorm()*v3.getNorm());
	float cos3=v3*v1/(v3.getNorm()*v1.getNorm());

	if(fabs(cos1)>1-eps||fabs(cos2)>1-eps||fabs(cos3)>1-eps)return true;
	else return false;
}

void HActiveContours::loadVertices( std::list<HPoint3<float>> &vertexList )
{

	if(vertices)delete vertices;
	vertices=new HVertexContainer(vertexList);

	//Cleaning all the auxiliary containers
	activeEdges.clear();
	frozenEdges.clear();
	triangles.clear();

	unvisitedCount=vertices->tree.size();
	return;

}

bool HActiveContours::growStep()
{	
	if(unvisitedCount==0&&activeEdges.size()==0)
	{
		unsigned int frozenBefore=frozenEdges.size();
		postStitch();
		return activeEdges.size()>0 || frozenEdges.size()<frozenBefore;
	}
	else if(activeEdges.size()>0)modelGrow();
	else modelInit();

	return true;
}

void HActiveContours::postStitch()
{
	if(frozenEdges.size()==0)return;

	for(std::list<HEdgeSeed>::iterator it=frozenEdges.begin(); it!=frozenEdges.end(); ++it )
	{
		HEdgeSeed e=*it;
		HPointSeed* pps1=getPropagatedVertex(e,true);
		if(!pps1)continue;

		bool isStitched=false;

		for(std::list<HEdgeSeed>::iterator ite=frozenEdges.begin(); ite!=frozenEdges.end(); ++ite )
		{
			HEdgeSeed ee=*ite;

			if(e==ee)continue;
			HPointSeed* pps2=getPropagatedVertex(ee,true);
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
					HPoint3<float> a=e.p2->p-e.p1->p;
					HPoint3<float> b=ee.p2->p-ee.p1->p;
					HPoint3<float> norm=getNormalVector(a,b);

					HPoint3<float> n1=getNormalVector(a,norm);
					float signCor1=n1*e.propagationVector>0?1.0f:-1.0f;
					n1=n1*signCor1;

					HPoint3<float> n2=getNormalVector(b,norm);
					float signCor2=n2*ee.propagationVector>0?1.0f:-1.0f;
					n2=n2*signCor2;

					float sgn1=a*n2;
					float sgn2=b*n1;
					float cosab=a*b/(a.getNorm()*b.getNorm());

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
							addActiveEdge(ne);

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


}