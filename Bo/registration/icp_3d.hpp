
/******************************************************************************

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
#include <boost/assert.hpp>

#include "bo/core/vector.hpp"
#include "bo/core/kdtree.hpp"
#include "bo/math/extended_math.hpp"
#include "bo/math/transformation_3d.hpp"
#include "bo/math/blas_extensions.hpp"
#include "bo/math/blas_conversions.hpp"

/* Usage example:

    std::vector<Vertex> std_source();
    ... // Initialize the source point cloud here.

    boost::shared_ptr<std::vector<Vertex> > std_target_ptr = boost::shared_ptr<std::vector<Vertex> >
    (new std::vector<Vertex>());
    ... // Initialize the target point cloud here.

    // Create an ICP registrator.
    bo::registration::ICP3D<float> icp(std_source, std_target_ptr,
    bo::distances::euclidean_distance<float, 3>);

    // Initialize auxiliary variables.
    // Iteration counter.
    std::size_t iteration = 0;
    // Norm of the transformation difference.
    float epsilon = std::numeric_limits<float>::max();
    // Distance between the source and target point clouds.
    float distance = std::numeric_limits<float>::max();
    // Transformation matrices.
    bo::math::matrix<float> m1;
    bo::math::matrix<float> m2;

    // Iteratively register two point clouds.
    while (iteration < max_allowed_iterations && epsilon >= min_transformation_epsilon)
    {
        // Cache the transformation from the previous iteration.
        m1 = m2;

        // Perform a new ICP registration step and update the current distance
        // between the source and target.
        distance = icp.next();

        // Update the current transformation.
        m2 = icp.current_transformation().matrix();

        // Update the norm of the transformation difference.
        if (iteration > 0)
        {
            bo::math::matrix<float> mdif = m2 - m1;
            epsilon = bo::math::l1_norm(mdif);
        }

        ++iteration;
    }

    // The final transformation.
    bo::Transformation3D<float> final_transformation = icp.current_transformation();

    // The corresponding final transformation matrix.
    bo::math::matrix<float> final_matrix = final_transformation.matrix();
*/

namespace bo {
namespace registration {

template <typename RealType>
class ICP3D: public boost::noncopyable
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> PointCloud;
    typedef boost::shared_ptr<PointCloud> PointCloudPtr;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
    typedef math::Transformation3D<RealType> Transformation;

public:
    ICP3D(const PointCloud& source, PointCloudPtr target, Metric dist_fun,
          bool is_preprocess = true);

    RealType next();
    Transformation current_transformation() const;
    const PointCloud& current_cloud() const;
    const PointCloud& current_correspondence() const;

private:

    typedef bo::KDTree<3, Point3D, std::pointer_to_binary_function
        <const Point3D&, std::size_t, RealType> > Point3DTree;
    static RealType point_bac(const Point3D& p, std::size_t k);

    void overlay_();
    Point3D centroid_(PointCloud* cloud) const;
    // Attention: the elements cloud2[i] must correspond to the elements cloud1[i].  
    math::matrix<RealType> cross_covariance_(PointCloud* cloud1, PointCloud* cloud2,
        const Point3D& centroid1, const Point3D& centroid2) const;
    // Attention: the elements cloud2[i] must correspond to the elements cloud1[i].  
    RealType distance_(PointCloud* cloud1, PointCloud* cloud2) const;
    void update_current_transform_and_cloud_(const Transformation& m);

private:
    PointCloudPtr target_cloud_;
    PointCloud current_cloud_;
    PointCloud corresp_cloud_;
    Transformation current_trans_;
    Point3DTree tree_;   
    Metric dist_fun_;
};


// Includes computing a kd-tree for the target point cloud.
template <typename RealType>
ICP3D<RealType>::ICP3D(const PointCloud& source, PointCloudPtr target,
    Metric dist_fun, bool is_preprocess):
    current_cloud_(source), target_cloud_(target), dist_fun_(dist_fun), 
    current_trans_(),
    tree_(target->begin(), target->end(), std::ptr_fun(ICP3D<RealType>::point_bac))
{    
    if (is_preprocess)
        overlay_();
}

template <typename RealType>
RealType ICP3D<RealType>::next()
{
    // Calculate the mass center of the current point cloud.
    Point3D current_centroid = centroid_(&current_cloud_);

    // Update the correspondence.
    corresp_cloud_.clear();
    for (typename PointCloud::const_iterator it = current_cloud_.begin(); it != current_cloud_.end(); ++it)
    {
        std::pair<typename Point3DTree::const_iterator, RealType> closest = tree_.find_nearest(*it);
        corresp_cloud_.push_back(*closest.first);
    }

    // Calculate the mass center of the corresponding points.
    Point3D corresp_centroid = centroid_(&corresp_cloud_);
    
    // Calculate the cross covariance for the current points and the target ones.
    math::matrix<RealType> Spx = cross_covariance_(&current_cloud_, &corresp_cloud_, current_centroid,
        corresp_centroid);

    math::matrix<RealType> SpxT = trans(Spx);

    // Create the asymmetrical matrix.
    math::matrix<RealType> Apx = Spx - SpxT;

    // Calculate the matrix trace.
    RealType traceSpx = Spx(0, 0) + Spx(1, 1) + Spx(2, 2);

    math::matrix<RealType> Bpx = Spx + SpxT -
        traceSpx * math::identity_matrix<RealType>(3);

    // Create the 4x4 matrix.
    math::matrix<RealType> Qpx(4, 4);
    
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
    math::subrange(Qpx, 1, 4, 1, 4) = math::subrange(Bpx, 0, 3, 0, 3);

    math::eigen_symmetric(Qpx);

    // Quaternion that defines the optimal rotation.
//    Vector<RealType, 4> quaternion = math::to_bo_vector(math::column(Qpx, 3));
    math::bounded_vector<RealType, 4> col(math::column(Qpx, 3));
    Vector<RealType, 4> quaternion = math::to_bo_vector(col);

    // Optimal translation. 
    Point3D translation = corresp_centroid - Transformation(quaternion) * current_centroid;

    // Create the optimal transformation.
    Transformation optimal_trans(quaternion, translation);

    // Update the current point cloud and the transformation.
    update_current_transform_and_cloud_(optimal_trans);

    return distance_(&current_cloud_, &corresp_cloud_);
}

template <typename RealType>
typename ICP3D<RealType>::Transformation ICP3D<RealType>::current_transformation() const
{
    return current_trans_;
}

template <typename RealType>
const typename ICP3D<RealType>::PointCloud& ICP3D<RealType>::current_cloud() const
{
    return current_cloud_;
}

template <typename RealType>
const typename ICP3D<RealType>::PointCloud& ICP3D<RealType>::current_correspondence() const
{
    return corresp_cloud_;
}

template <typename RealType>
void ICP3D<RealType>::overlay_()
{
    // Cloud shift using centroids.

    Point3D translation = centroid_(target_cloud_.get()) - centroid_(&current_cloud_);

    Transformation shift_transform(translation);
    
    // Update the current point cloud and the transformation.
    update_current_transform_and_cloud_(shift_transform);
}

template <typename RealType>
typename ICP3D<RealType>::Point3D ICP3D<RealType>::centroid_(PointCloud* cloud) const
{
    BOOST_ASSERT(cloud->size() > 0);
    Point3D mass_center = math::mean(*cloud);

    return mass_center;
}

template <typename RealType>
math::matrix<RealType> ICP3D<RealType>::cross_covariance_(PointCloud* cloud1, PointCloud* cloud2,
    const Point3D& centroid1, const Point3D& centroid2) const
{
    const std::size_t n = cloud1->size();

    BOOST_ASSERT(cloud2->size() == n && n > 0);

    math::matrix<RealType> m = math::zero_matrix<RealType>(3, 3);
    math::matrix<RealType> p(3, 1);
    math::matrix<RealType> x(1, 3);

    for (std::size_t i = 0; i < n; ++i)
    {
        Point3D a = cloud1->at(i);
        Point3D b = cloud2->at(i);

        p(0, 0) = a[0] - centroid1[0];
        p(1, 0) = a[1] - centroid1[1];
        p(2, 0) = a[2] - centroid1[2];

        x(0, 0) = b[0] - centroid2[0];
        x(0, 1) = b[1] - centroid2[1];
        x(0, 2) = b[2] - centroid2[2];

        // Accumulate the elements.
        // Warning: type overflow is possible here!
        m = m + math::prod(p, x);
    }   

    m = m / n;

    return m;
}

template <typename RealType>
RealType ICP3D<RealType>::distance_(PointCloud* cloud1, PointCloud* cloud2) const
{   
    const std::size_t n = cloud1->size();

    BOOST_ASSERT(cloud2->size() == n);

    RealType sum(0);

    for (std::size_t i = 0; i < n; ++i)
    {
        Point3D a = cloud1->at(i);
        Point3D b = cloud2->at(i);
        sum += dist_fun_(a, b);
    }

    return sum;
}

// Point3D brackets accessor.
template <typename RealType>
inline RealType ICP3D<RealType>::point_bac(const Point3D& p, std::size_t k)
{ 
    return p[k]; 
}

template <typename RealType>
void ICP3D<RealType>::update_current_transform_and_cloud_(const Transformation& m)
{
    // Update the current transformation.
    current_trans_ = m * current_trans_;

    // Update the current point cloud.
    for (typename PointCloud::iterator it = current_cloud_.begin(); it != current_cloud_.end(); ++it)
    {
        *it = m * (*it);
    }
}

} // namespace registration
} // namespace bo

#endif // ICP_3D_HPP_507CF525_CC8A_4499_80D2_E8C8644603F5_
