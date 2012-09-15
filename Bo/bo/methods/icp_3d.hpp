
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
    typedef boost::array<boost::array<RealType, 4> > Matrix4x4;

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
    Matrix4x4 cross_covariance_(PointCloud* cloud1, PointCloud* cloud2,
        Point3D centroid1, Point3D centroid2, Correspondences* corresp) const;

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
    current_cloud_(source), target_cloud(target), dist_fun_(dist_fun)
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
    RealType error(0);

    // TODO: implement ICP iteration here.

    return error;
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
    // TODO: return centroid for a given point cloud.
}

template <typename RealType>
typename ICP3D<RealType>::Matrix4x4 ICP3D<RealType>::cross_covariance_(PointCloud* cloud1,
    PointCloud* cloud2, Point3D centroid1, Point3D centroid2, Correspondences* corresp) const
{
    // TODO: return a 4x4 cross-covariance matrix for given point clouds.
}

} // namespace methods
} // namespace bo

#endif // ICP_3D_HPP_507CF525_CC8A_4499_80D2_E8C8644603F5_
