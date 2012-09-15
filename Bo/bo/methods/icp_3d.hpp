
/******************************************************************************

  icp_3d.hpp, v 0.0.1 2012.09.15

  Point-to-point implementation of the ICP registration algorithm for 3D point
  clouds.

  Copyright (c) 2012
  Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef ICP_3D_HPP_507CF525_CC8A_4499_80D2_E8C8644603F5_
#define ICP_3D_HPP_507CF525_CC8A_4499_80D2_E8C8644603F5_

#include <vector>
#include <utility>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>

#include "bo/vector.hpp"
#include "bo/blas/blas.hpp"

namespace bo {
namespace methods {

namespace detail {

} // namespace detail

template <typename RealType>
class ICP3D: public boost::noncopyable
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> PointCloud;
    typedef boost::shared_ptr<PointCloud> PointCloudPtr;
    typedef std::pair<std::size_t, std::size_t> Correspondence;
    typedef std::vector<Correspondence> Correspondences;
    typedef boost::shared_ptr<Correspondences> CorrespondencesPtr;
    typedef boost::function<RealType (Point3D, Point3D)> DistanceFunction;

public:
    ICP3D(const PointCloud& source, PointCloudPtr target, DistanceFunction dist_fun,
          bool is_preprocess = true);

    RealType next();
    Transformation3D current_transformation() const;
    PointCloudPtr current_cloud() const;
    CorrespondencesPtr current_correspondence() const;

private:
    void overlay_();
    Point3D centroid_(PointCloud* cloud) const;
    matrix<RealType> cross_covariance_(PointCloud* cloud1, PointCloud* cloud2,
        const Point3D& centroid1, const Point3D& centroid2, Correspondences* corresp) const;
    RealType distance_(PointCloud* cloud1, PointCloud* cloud2, Correspondences* corresp) const;

private:
    PointCloudPtr target_cloud_;
    Point3D target_centroid_;
    DistanceFunction dist_fun_;
    PointCloud current_cloud_;
    CorrespondencesPtr current_corresp_;
    Transformation3D current_trans_;
};

template <typename RealType>
ICP3D<RealType>::ICP3D(const PointCloud& source, PointCloudPtr target,
    DistanceFunction dist_fun, bool is_preprocess):
    current_cloud_(source), target_cloud(target), dist_fun_(dist_fun), 
    current_trans_()
{
    // Cache centroid for the target point cloud.
    target_centroid_ = centroid_(target);

    // Precompute a kd-tree for the target point cloud.

    // TODO: prepare for ICP iterations.

    if (is_preprocess)
        overlay_();
}

template <typename RealType>
RealType ICP3D<RealType>::next()
{
    Point3D current_centroid = centroid_(current_cloud_);

    matrix<RealType> Spx = cross_covariance_(current_cloud_, target_cloud_, current_centroid,
        target_centroid_, current_corresp_);

    matrix<RealType> SpxT = trans(Spx);

    // Create the asymmetrical matrix.
    matrix<RealType> Apx = Spx - SpxT;

    // Calculate the matrix trace.
    RealType traceSpx = Spx(0, 0) + Spx(1, 1) + Spx(2, 2);

    matrix<RealType> Bpx = Spx + SpxT + traceSpx * identity_matrix<RealType>(3);

    // Create the 4x4 matrix.
    matrix<RealType> Qpx(4, 4);
    
    // Fill in the matrix.
    // Block 1.
    Qpx(0, 0) = traceSpx;
    // Block 2.
    Qpx(0, 1) = Apx(1, 2); 
    Qpx(0, 2) = Apx(2, 0);
    Qpx(0, 3) = Apx(0, 1);
    // Block 3.
    Qpx(1, 0) = Apx(1, 2); 
    Qpx(2, 0) = Apx(2, 0);
    Qpx(3, 0) = Apx(0, 1);
    // Block 4.
    Qpx(1, 1) = Bpx(0, 0);
    Qpx(1, 2) = Bpx(0, 1);
    Qpx(1, 3) = Bpx(0, 2);
    Qpx(2, 1) = Bpx(1, 0);
    Qpx(2, 2) = Bpx(1, 1);
    Qpx(2, 3) = Bpx(1, 2);
    Qpx(3, 1) = Bpx(2, 0);
    Qpx(3, 2) = Bpx(2, 1);
    Qpx(3, 3) = Bpx(2, 2);

    eigen_analysis(Qpx);

    // Quaternion that defines the optimal rotation.
    Vector<RealType, 4> quaternion;
    quaternion[0] = Qpx(3, 0);
    quaternion[1] = Qpx(3, 1);
    quaternion[2] = Qpx(3, 2);
    quaternion[3] = Qpx(3, 3);
    Vector<RealType, 3> t(0);

    // Optimal translation. 
    Point3D translation = target_centroid_ - Transformation3D(quaternion, t) * current_centroid;

    // Create the optimal transformation.
    Transformation3D optimal_trans(quaternion, translation);

    // Update the current transformation.
    current_trans_ = optimal_trans * current_trans_;

    // Update the current point cloud.
    PointCloud::iterator it = current_cloud_.begin();
    while (it != current_cloud_.end())
    {
        *it = current_trans_ * (*it);
        ++it;
    }

    return distance_(current_cloud_, target_cloud_, current_corresp_);
}

template <typename RealType>
Transformation3D ICP3D<RealType>::current_transformation() const
{
    return current_transformation_;
}

template <typename RealType>
typename ICP3D<RealType>::PointCloudPtr ICP3D<RealType>::current_cloud() const
{
    return current_cloud_;
}

template <typename RealType>
typename ICP3D<RealType>::CorrespondencesPtr ICP3D<RealType>::current_correspondence() const
{
    return current_corresp_;
}

template <typename RealType>
void ICP3D<RealType>::overlay_()
{
    // TODO: implement current cloud shift using centroids.
}

template <typename RealType>
typename ICP3D<RealType>::Point3D ICP3D<RealType>::centroid_(PointCloud* cloud) const
{
    Point3D mass_center(0);
    
    PointCloud::const_iterator it = cloud->begin();
    
    while (it != cloud->end())
    {
        mass_center += *it;    
        ++it;
    }

    BOOST_ASSERT(cloud->size() != 0);

    mass_center /= cloud->size();

    return mass_center;
}

template <typename RealType>
matrix<RealType> ICP3D<RealType>::cross_covariance_(PointCloud* cloud1, PointCloud* cloud2,
    const Point3D& centroid1, const Point3D& centroid2, Correspondences* corresp) const
{
    matrix<RealType> m(3, 3) = zero_matrix<RealType>(3, 3);
    matrix<RealType> p(3, 1);
    matrix<RealType> x(1, 3);   

    Correspondences::const_iterator it = corresp->begin();

    // Calculate cross-covariance.
    while (it != corresp->end())
    {       
        Point3D a = cloud1->at(it->first);
        Point3D b = cloud2->at(it->second);

        p(0, 0) = a[0] - centroid1[0];
        p(1, 0) = a[1] - centroid1[1];
        p(2, 0) = a[2] - centroid1[2];

        x(0, 0) = b[0] - centroid2[0];
        x(0, 1) = b[1] - centroid2[1];
        x(0, 2) = b[2] - centroid2[2];

        m = m + prod(p, x);
        
        ++it;
    }

    BOOST_ASSERT(corresp->size() != 0);

    m = m / corresp->size();

    return m;
}

template <typename RealType>
RealType bo::methods::ICP3D<RealType>::distance_( PointCloud* cloud1, PointCloud* cloud2, 
    Correspondences* corresp ) const
{
    RealType sum(0);

    Correspondences::const_iterator it = corresp->begin();

    // Calculate the integral distance.
    while (it != corresp->end())
    {       
        Point3D a = cloud1->at(it->first);
        Point3D b = cloud2->at(it->second);

        sum += dist_fun_(a, b);

        ++it;
    }

    return sum;
}

} // namespace methods
} // namespace bo

#endif // ICP_3D_HPP_507CF525_CC8A_4499_80D2_E8C8644603F5_
