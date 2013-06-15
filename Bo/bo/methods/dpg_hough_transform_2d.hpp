
/******************************************************************************

  Dual-point generalized Hough Transform for 2D object recognition.

  Copyright (c) 2013
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

#ifndef DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419_
#define DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419_

#include <utility>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <algorithm>
#include <complex>

#include <boost/noncopyable.hpp>
#include <boost/assert.hpp>
#include <boost/math/constants/constants.hpp>

#include "bo/vector.hpp"
#include "bo/blas/blas.hpp"
#include "bo/raw_image_2d.hpp"
#include "bo/topology.hpp"
#include "bo/blas/blas.hpp"
#include "bo/methods/convex_hull_3d.hpp"

#include <iostream>

namespace bo {
namespace methods {
namespace recognition {

namespace detail{

// A 4-dimensional hyperplane defined by a point located
// on this plane and the plane's normal vector.
template <typename RealType>
class Hyperplane4D
{
public:
    typedef Vector<RealType, 4> Point4D;

    Hyperplane4D(const Point4D &point, const Point4D &normal)
        : point_(point), kEpsilon(0.0001)
    {
        Point4D n = normal / normal.euclidean_norm();
        normal_ = n;

        // Decomposition matrix (inverse of the base one).
        m_ = blas::matrix<RealType>(4, 4);
        m_(0, 0) =  n[3]; m_(0, 1) =  n[2]; m_(0, 2) = -n[1]; m_(0, 3) = -n[0];
        m_(1, 0) = -n[2]; m_(1, 1) =  n[3]; m_(1, 2) =  n[0]; m_(1, 3) = -n[1];
        m_(2, 0) =  n[1]; m_(2, 1) = -n[0]; m_(2, 2) =  n[3]; m_(2, 3) = -n[2];
        m_(3, 0) =  n[0]; m_(3, 1) =  n[1]; m_(3, 2) =  n[2]; m_(3, 3) =  n[3];
    }

    inline const Point4D& point() const
    {
        return point_;
    }

    inline const Point4D& normal() const
    {
        return normal_;
    }

    // Computes intersection of a line defined by two 4-dimensional points
    // with the hyperplane. Considers the intersection as a linear coordinate 't'
    // such that the intersection point is calculated as P = q1 + t * (q2 - q1).
    // If the line intersects with the hyperplane, returns true and updates the
    // function parameter 't'. Otherwise returns false.
    bool intersect(const Point4D& q1, const Point4D& q2, RealType& t)
    {
        Point4D q1p = point_ - q1;
        Point4D e = q2 - q1;
        RealType dot1 = q1p * normal_;

        RealType e_norm = e.euclidean_norm();
        // This case is needed for continuity, when two very close points
        // both belong to the hyperplane.
        if (e_norm < kEpsilon)
        {
            if (dot1 < kEpsilon)
            {
                t = RealType(0);
                return true;
            }
            else
            {
                return false;
            }
        }

        RealType dot2 = e * normal_;

        if (std::abs(dot2 / e_norm) > kEpsilon)
        {
            // If the line is not parallel to the hyperplane.
            t = dot1 / dot2;
            return true;
        }
        else
        {
            return false;
        }
    }

    // Computes the decomposition of the point in the hyperplane's
    // coordinate system defined by the decomposition matrix.
    Point4D decompose(const Point4D& p)
    {
        Point4D q = p - point_;

        blas::matrix<RealType> v(4, 1);
        v(0, 0) = q[0];
        v(1, 0) = q[1];
        v(2, 0) = q[2];
        v(3, 0) = q[3];

        // Compute decomposition in the projection basis.
        v = blas::prod(m_, v);

        return Point4D(v(0, 0), v(1, 0), v(2, 0), v(3, 0));
    }

private:

    Point4D point_;
    Point4D normal_;

    // Hyperplane's decomposition matrix.
    blas::matrix<RealType> m_;

    const RealType kEpsilon;
};

template <typename RealType>
class Space
{
public:

    typedef Vector<RealType, 4> Point4D;
    typedef std::vector<Point4D> Points4D;
    typedef Vector<std::size_t, 4> Size4D;
    typedef std::pair<Point4D, Point4D> Box4D;
    typedef bo::topology::OrthotopeTopology<RealType, 4> Geometry;
    typedef detail::Hyperplane4D<RealType> Hyperplane;

    // A line in 4D is defined by any point located
    // on this line and its directional vector.
    struct Line4D
    {
        Line4D()
            : point(Point4D(0)), direction(Point4D(0))
        { }

        Line4D(const Point4D &_point, const Point4D &_direction)
            : point(_point), direction(_direction)
        { }

        Point4D point;
        Point4D direction;
    };



    typedef Vector<RealType, 2> SegmentCoordinates;
    typedef std::vector<Space> Spaces;

    // A segment in 4D is modelled by a 4D line and the
    // coordinates of two points on this line
    // relatively to the line's direction vector.
    typedef std::pair<Line4D, SegmentCoordinates> Segment4D;

    Space(const Box4D &box = Box4D(Point4D(0, 0, 0, 0), Point4D(0, 0, 0, 0)),
          const Size4D &divisions_per_dimension = Size4D(2, 2, 2, 2),
          std::size_t max_resolution_level = 1,
          std::size_t cell_resolution_increment = 1,
          std::size_t resolution_level = 0):
        box_(box),  divisions_per_dimension_(divisions_per_dimension),
        max_resolution_level_(max_resolution_level), resolution_level_(resolution_level),
        votes_(0)
    {
        BOOST_ASSERT(resolution_level_ <= max_resolution_level_);

        // Adjust the cell resolution increment.
        if (resolution_level_ + cell_resolution_increment > max_resolution_level_)
        {
            cell_resolution_increment_ = max_resolution_level_ - resolution_level_;
        }
        else
        {
            cell_resolution_increment_ =  cell_resolution_increment;
        }

        // Compute the size of cells used for vote calculation.
        for (std::size_t i = 0; i < 4; ++i)
        {
            std::size_t cells_in_dimension = static_cast<std::size_t>(
                std::pow(RealType(divisions_per_dimension_[i]),
                         RealType(max_resolution_level_) - RealType(resolution_level_)));

            cell_size_[i] = (box_.second[i] - box_.first[i]) /  cells_in_dimension;
        }

        min_cell_volume_ = std::pow((cell_size_[0] + cell_size_[1] + cell_size_[2] + cell_size_[3]) / 4, 3);

        // Compute the number of probabilistic elements used for the subdivision policy.
        prob_element_count_ = 1;
        for (std::size_t i = 0; i < 4; ++i)
        {
            std::size_t prob_elements_in_dimension = static_cast<std::size_t>(
                std::pow(RealType(divisions_per_dimension_[i]), RealType(cell_resolution_increment_)));

            prob_element_count_ *= prob_elements_in_dimension;
        }
    }

    // Subdivides the space.
    void subdivide()
    {
        subspaces_.clear();

        std::size_t subdivision_count_ = divisions_per_dimension_[0] * divisions_per_dimension_[1] *
                                         divisions_per_dimension_[2] * divisions_per_dimension_[3];

        subspaces_.reserve(subdivision_count_);

        Point4D sizes4d = box_.second - box_.first;
        Point4D steps4d = Point4D(sizes4d[0] / RealType(divisions_per_dimension_[0]),
                                  sizes4d[1] / RealType(divisions_per_dimension_[1]),
                                  sizes4d[2] / RealType(divisions_per_dimension_[2]),
                                  sizes4d[3] / RealType(divisions_per_dimension_[3]));

        for (std::size_t d0 = 0; d0 < divisions_per_dimension_[0]; ++d0)
            for (std::size_t d1 = 0; d1 < divisions_per_dimension_[1]; ++d1)
                for (std::size_t d2 = 0; d2 < divisions_per_dimension_[2]; ++d2)
                    for (std::size_t d3 = 0; d3 < divisions_per_dimension_[3]; ++d3)
                    {
                        // Create a subspace.
                        Point4D translation(d0 * steps4d[0], d1 * steps4d[1], d2 * steps4d[2], d3 * steps4d[3]);
                        Box4D b(box_.first + translation, box_.first + translation + steps4d);
                        Space s(b, divisions_per_dimension_, max_resolution_level_, cell_resolution_increment_,
                                resolution_level_ + 1);

                        subspaces_.push_back(s);
                    }

    }

    inline Spaces& get_subspaces()
    {
        return subspaces_;
    }

    inline const Spaces& get_subspaces() const
    {
         return subspaces_;
    }

    inline const Box4D& get_bounding_box() const
    {
        return box_;
    }

    inline RealType get_votes() const
    {
        return votes_;
    }

    inline std::size_t get_prob_element_count() const
    {
        return prob_element_count_;
    }

    void reset_votes()
    {
        votes_ = 0;
    }

    inline std::size_t get_resolution_level() const
    {
        return resolution_level_;
    }

    inline std::size_t get_max_resolution_level() const
    {
        return max_resolution_level_;
    }

    inline Point4D get_mass_center() const
    {
        return (box_.first + box_.second) / 2;
    }

    void vote(Segment4D segment)
    {
        // The number of votes that the space receives in the result of the intersection
        // equals to the lenght of the resulting segment.
        // Compute the segment.
        Point4D s = segment.first.direction * (segment.second[1] - segment.second[0]);

        // Increase the votes.
        // Continuous. Attention: has no effect if the intersection is in one point.
        votes_ += s.euclidean_norm();
    }

    void vote_unit(const Segment4D &segment)
    {
        // Binary.
        if (segment.second[0] > -std::numeric_limits<RealType>::max() &&
            segment.second[0] < std::numeric_limits<RealType>::max())
        {
            votes_ += 1;
        }
    }

    // Increases the number of votes on the number of intersection of the given segment
    // with the grid defined by the space cell size.
    void vote_descrete(const Segment4D &segment)
    {
        // We will accumulate the coordinates of intersections of the given segment
        // with the grid of the given cell size.
        std::set<RealType> intersections;

        // Find coordinates of the intersection of the line with the space.
        Segment4D base_seg = intersect(segment.first);


        RealType t1_base = base_seg.second[0];
        RealType t2_base = base_seg.second[1];
        // Order.
        order(t1_base, t2_base);

        RealType t1 = segment.second[0];
        RealType t2 = segment.second[1];
        // Order.
        order(t1, t2);

        // Segment direction.
        Point4D v = segment.first.direction;

        // For each dimension define the coordinate step beetween two cells (dimension levels).
        for (std::size_t d = 0; d < 4; ++d)
        {
            if (v[d] != 0)
            // If it is zero, the segment does not intersect the levels of this dimension and
            // we can skip further analysis.
            {
                RealType dt = std::abs(cell_size_[d] / v[d]);

                // Accumulate all intersections of this dimension situated between the begin
                // and the end of the segment. TODO: optimize it!
                for (RealType t = t1_base + dt; t < t2_base; t += dt)
                {
                    if (t > t1 && t < t2)
                        intersections.insert(t);
                }
            }
        }

        votes_ += intersections.size() + 1;
    }

    // Increases the votes on the number equal to the taxicab norm of the
    // corresponding inner segment relatively to the space grid.
    void vote_taxicab(const Segment4D &segment)
    {
        // If the segment is not zero.
        if (segment.second[0] == segment.second[1])
            return;

        std::size_t taxicab_dst = 0;

        // Local coordinate origin.
        Point4D bot = box_.first;
        Point4D top = box_.second;

        // The end points of the segment.
        Point4D p1 = segment.first.point + segment.first.direction * segment.second[0];
        Point4D p2 = segment.first.point + segment.first.direction * segment.second[1];

        // Segment direction.
        Point4D v = segment.first.direction;

        // For each dimension define the coordinate step beetween two cells (dimension levels).
        for (std::size_t d = 0; d < 4; ++d)
        {
            if (v[d] != 0)
            // If it is zero, the segment does not intersect the levels of this dimension and
            // we can skip further analysis.
            {
                RealType proj1 = p1[d];
                RealType proj2 = p2[d];
                order(proj1, proj2);

                order(bot[d], top[d]);

                cut(proj1, bot[d], top[d]);
                cut(proj2, bot[d], top[d]);

                BOOST_ASSERT(proj1 >= bot[d] && proj2 >= bot[d]);

                RealType block1 = std::floor((proj1 - bot[d]) / cell_size_[d]);
                RealType block2 =  std::ceil((proj2 - bot[d]) / cell_size_[d]);

                taxicab_dst += static_cast<std::size_t>(block2 - block1);
            }
        }

        votes_ += taxicab_dst;
    }

    // Increases the votes on the number equal to the maximum (uniform) norm of the
    // inner segment relatively to the space grid.
    void vote_maxnorm(const Segment4D &segment)
    {
        // If the segment is not zero.
        if (segment.second[0] == segment.second[1])
            return;

        std::size_t maximum_norm = 0;

        // Local coordinate origin.
        Point4D bot = box_.first;
        Point4D top = box_.second;

        // The end points of the segment.
        Point4D p1 = segment.first.point + segment.first.direction * segment.second[0];
        Point4D p2 = segment.first.point + segment.first.direction * segment.second[1];

        // Segment direction.
        Point4D v = segment.first.direction;

        // For each dimension define the coordinate step beetween two cells (dimension levels).
        for (std::size_t d = 0; d < 4; ++d)
        {
            if (v[d] != 0)
            // If it is zero, the segment does not intersect the levels of this dimension and
            // we can skip further analysis.
            {
                RealType proj1 = p1[d];
                RealType proj2 = p2[d];
                order(proj1, proj2);

                order(bot[d], top[d]);

                cut(proj1, bot[d], top[d]);
                cut(proj2, bot[d], top[d]);

                BOOST_ASSERT(proj1 >= bot[d] && proj2 >= bot[d]);

                RealType block1 = std::floor((proj1 - bot[d]) / cell_size_[d]);
                RealType block2 =  std::ceil((proj2 - bot[d]) / cell_size_[d]);

                std::size_t dst = static_cast<std::size_t>(block2 - block1);

                if (maximum_norm < dst)
                    maximum_norm = dst;
            }
        }

        votes_ += maximum_norm;
    }

    // Returns the result of intersection of the given segment with the space.
    Segment4D intersect(Segment4D segment)
    {
        // Cut the segment in each of four dimensions.
        for (std::size_t d = 0; d < 4; ++d)
            cut_segment(segment, d, box_.first[d], box_.second[d]);

        return segment;
    }

    // Returns the result of intersection (segment) of the given 4D line with the space.
    Segment4D intersect(Line4D line)
    {
        // Create the "infinite" segment coordinates that models the whole line.
        RealType mmax = std::numeric_limits<RealType>::max();
        SegmentCoordinates coord(-mmax, mmax);

        // Create the "infinite" segment.
        Segment4D infinite_seg(line, coord);

        return intersect(infinite_seg);
    }

    // Returns the points that define the polyhedron of the hyperplane-hyperrectangle
    // intersection.
    Points4D intersect(Hyperplane plane)
    {
        Points4D vertices;

        Point4D d = box_.second - box_.first;

        typename Geometry::Edges edges = Geometry::edges();

        for (typename Geometry::Edges::const_iterator it = edges.begin();
             it != edges.end(); ++it)
        {
            typename Geometry::Point e1 = it->first;
            typename Geometry::Point e2 = it->second;

            // Calculate the adjacent vertices of the box edges.
            Point4D q1 = box_.first + Point4D(e1[0] * d[0], e1[1] * d[1],
                                              e1[2] * d[2], e1[3] * d[3]);
            Point4D q2 = box_.first + Point4D(e2[0] * d[0], e2[1] * d[1],
                                              e2[2] * d[2], e2[3] * d[3]);

            // Calculate intersection of the edge with the plane.
            RealType t;
            if (plane.intersect(q1, q2, t))
            {
                // Add the intersection point if it belongs to the edge.
                if (t >= 0 && t <= 1)
                {
                    vertices.push_back(q1 + t * (q2 - q1));
                }
            }
        }

        return vertices;
    }

    RealType intersection_volume(Hyperplane plane)
    {
        RealType kEpsilon(0.001);

        // The points of intersection with the plane in 4D.
        Points4D vertices4 = intersect(plane);

        // Container for projections.
        typedef methods::surfaces::IncrementalConvexHull3D<RealType> Hull;
        typename Hull::Points3D vertices3;

        // Project the vertices into 3D.
        for (typename Points4D::const_iterator it = vertices4.begin();
             it != vertices4.end(); ++it)
        {
            Point4D p = plane.decompose(*it);

            // The vertex must lie in 3D!
            BOOST_ASSERT(std::abs(p[3]) < kEpsilon);

            // Projection.
            typename Hull::Point3D pr(p[0], p[1], p[2]);

            vertices3.push_back(pr);
        }

        if(vertices3.size() > 3)
        {
            // Compute the convex hull and its volume.
            Hull hull3d(vertices3);
            return hull3d.get_volume();
        }
        else
        {
            return 0;
        }
    }

    // The number of votes that the space receives is equal to the volume
    // of the resulting hyperrectangle-hyperplane intersection.
    void vote(Hyperplane plane)
    {
        RealType volume = intersection_volume(plane);

        votes_ += volume;
    }

    // The number of votes that the space receives is equal to the number
    // of the cells that the hyperrectangle-hyperplane intersection "covers".
    void vote_descrete(Hyperplane plane)
    {
        RealType volume = intersection_volume(plane);

        RealType vote = volume / min_cell_volume_;

        votes_ += std::floor(vote) + 1;
    }

    void vote_unit(Hyperplane plane)
    {
        if (intersect(plane).size() > 0)
            votes_ += 1;
    }

private:

    Spaces subspaces_;
    Box4D box_;
    Size4D divisions_per_dimension_;
    std::size_t max_resolution_level_;
    std::size_t cell_resolution_increment_;
    std::size_t resolution_level_;
    RealType votes_;
    Point4D cell_size_;
    RealType min_cell_volume_;
    std::size_t prob_element_count_;

    // Cuts the segment in the given dimension according to two given levels.
    void cut_segment(Segment4D &seg, std::size_t dimension, RealType level1, RealType level2)
    {
        // Sort the levels: level1 <= level2.
        order(level1, level2);

        // The origin point and the direction.
        Point4D p = seg.first.point;
        Point4D v = seg.first.direction;

        // If the segment is parallel to the levels of this dimension.
        if (v[dimension] == 0)
        {
             // If the origin is outside the levels, return infinite point.
            if (p[dimension] > level2 || p[dimension] < level1)
                seg.second[0] = seg.second[1] = std::numeric_limits<RealType>::max();

            // Otherwise return the same segment.
            return;
        }

        // Define the coordinates of the line instersection with the dimension levels.
        RealType t1 = (level1 - p[dimension]) / v[dimension];
        RealType t2 = (level2 - p[dimension]) / v[dimension];

        // Sort the coordinates: t1 <= t2.
        order(t1, t2);

        // If the segments is outside the levels (less), move it in the negative
        // infinite point.
        if (seg.second[0] < t1 && seg.second[1] < t1)
        {
            seg.second[0] = seg.second[1] = -std::numeric_limits<RealType>::max();
            return;
        }

        // If the segments is outside the levels (less), move it in the positive
        // infinite point.
        if (seg.second[0] > t2 && seg.second[1] > t2)
        {
            seg.second[0] = seg.second[1] = std::numeric_limits<RealType>::max();
            return;
        }

        // Cut the segment.
        for (int i = 0; i < 2; ++i)
        {
            if (seg.second[i] < t1)
                seg.second[i] = t1;

            if (seg.second[i] > t2)
                seg.second[i] = t2;
        }

    }

    inline void order(RealType &a, RealType &b)
    {
        if (a > b) std::swap(a, b);
    }

    inline void cut(RealType &x, RealType a, RealType b)
    {
        if (x < a)
            x = a;
        if (x > b)
            x = b;
    }
};



template <typename RealType>
class SubdivisionPolicy
{
public:
    typedef Space<RealType> Space4D;

    // Returns the minimal number of spaces (at the current resolution level) from the beginning of the given
    // sorted collection such that occurrence of the subspaces (at the cell resolution level) with maximal
    // number of votes is not less then the given probability p.
    // Attention: the input collection of spaces must be sorted in descending order!
    static std::size_t probabilistic(const typename Space4D::Spaces &spaces, RealType p)
    {
        typename Space4D::Spaces subcollection1;
        typename Space4D::Spaces subcollection2 = spaces;

        std::size_t n = 0;

        // Find the minimal number of spaces that satisfy the probability constraint.
        while (n < spaces.size() &&
               p_max_value_in_subcollection(subcollection1, subcollection2) < p)
        {
            typename Space<RealType>::Spaces::iterator it = subcollection2.begin();

            subcollection1.push_back(*it);
            subcollection2.erase(it);
            n = subcollection1.size();
        }

        return n;
    }


    // Computes the probability of the event that the maximal value of subcollection1 is greater
    // than the maximal value from subcollection2.
    static RealType p_max_value_in_subcollection(const typename Space4D::Spaces &subcollection1,
                                                 const typename Space4D::Spaces &subcollection2)
    {
        // Compute the maximal vote in subcollection1.
        std::size_t max_votes = 0;
        for (typename Space4D::Spaces::const_iterator it = subcollection1.begin();
             it != subcollection1.end(); ++it)
        {
            std::size_t votes = std::size_t(it->get_votes());
            if (max_votes < votes)
            {
                max_votes = votes;
            }
        }

        RealType p = 0;

        // Compute the probability.
        for (std::size_t t = 1; t <= max_votes; ++t)
        {
            RealType f1 = F_joint(subcollection1, t);
            RealType f2 = F_joint(subcollection1, t - 1);
            RealType f3 = F_joint(subcollection2, t - 1);

            p += (f1 - f2) * f3;
        }

        return p;
    }

    // Joint distribution function for spaces.
    static RealType F_joint(const typename Space4D::Spaces &spaces, std::size_t x)
    {
        if (spaces.size() == 0) return 0;

        RealType f = 1;

        for (typename Space4D::Spaces::const_iterator it = spaces.begin();
             it != spaces.end(); ++it)
        {
            f *= F(std::size_t(it->get_votes()), it->get_prob_element_count(), x);
        }

        return f;
    }

    // Distribution function.
    static RealType F(std::size_t k, std::size_t n, std::size_t x)
    {
        // Some cases of small values (k, n) are implemented explicitly
        // avoiding the Gumbel-based approximation in order to increase
        // precision.
        if (k == 0)
        {
            return RealType(1);
        }
        if (n == 1)
        {
            return (x >= k) ? RealType(1) : RealType(0);
        }

        return std::exp(-std::exp((mu(k, n) - x) / beta(k, n)));
    }

    // Mean.
    static RealType E(std::size_t k, std::size_t n)
    {
        // Euler–Mascheroni constant.
        const RealType gamma = RealType(0.5772);

        return mu(k, n) + gamma * beta(k, n);
    }

    static RealType mu(std::size_t k, std::size_t n)
    {
        return RealType(1.16) * k * std::pow(RealType(n), RealType(-2) / 3) + 1;
    }

    static RealType beta(std::size_t k, std::size_t n)
    {
        return RealType(0.4) * k * std::pow(RealType(n), RealType(-4) / 5) + RealType(0.32);
    }
};



// Compares two spaces using an approximation of the Gumbold means.
template <typename RealType>
bool operator < (const Space<RealType> &s1, const Space<RealType> &s2)
{
    RealType mean1 = SubdivisionPolicy<RealType>::E(std::size_t(s1.get_votes()),
                                                    s1.get_prob_element_count());
    RealType mean2 = SubdivisionPolicy<RealType>::E(std::size_t(s2.get_votes()),
                                                    s2.get_prob_element_count());
    return
        (mean1 < mean2) ? true : false;
}



} // namespace detail

template <typename RealType>
class DualPointGHT: public boost::noncopyable
{
public:
    typedef DualPointGHT<RealType> this_type;
    typedef Vector<RealType, 2> Point2D;
    typedef std::vector<Point2D> Points2D;
    typedef std::pair<Point2D, Point2D> Reference;
    typedef std::pair<Point2D, Point2D> SearchArea;
    typedef std::pair<Reference, RealType> ReferenceVote;
    typedef std::vector<ReferenceVote> ReferenceVotes;
    // Feature is a model point and a tangent vector.
    typedef std::pair<Point2D, Point2D> Feature;
    typedef std::vector<Feature> Features; 
    typedef std::pair<RealType, RealType> ATableElement;
    typedef std::list<ATableElement> ATableRow;
    typedef std::vector<ATableRow> ATable;
    typedef detail::Space<RealType> Space4D;
    typedef detail::SubdivisionPolicy<RealType> SubPolicy;

    DualPointGHT(const Features &model_features, const Reference &model_reference, 
                 RealType tangent_accuracy = RealType(0.005)):
    model_reference_(model_reference), tangent_accuracy_(tangent_accuracy),
    pi_(boost::math::constants::pi<RealType>())
    {
        encode(model_features);

        model_base_ = (model_reference_.second - model_reference_.first).euclidean_norm();
    }

    void project_detected_lines(const Features &object_features, const Point2D &scaling_range,
                                bo::RawImage2D<RealType> &image1, bo::RawImage2D<RealType> &image2)
    {
        for (typename Features::const_iterator it = object_features.begin();
             it != object_features.end(); ++it)
        {
            // Reconstruct all the 4D lines from the alpha-table relatively to the current feature.
            for (std::size_t index = 0; index < atable_.size(); ++index)
            {
                // Probable angle between the directional vector and the tangent.
                RealType gamma = atable_gamma(index);

                // Find all corresponding directional vectors v2 defined by the (alpha, beta) angles of
                // the current table row.
                for (typename ATableRow::const_iterator abit = atable_.at(index).begin();
                     abit !=  atable_.at(index).end(); ++abit)
                {
                    typename Space4D::Line4D line4 = line4_from_feature_and_atable_element(*it, gamma, *abit);

                    // Create two segments on this line that correspond to the given scaling range.
                    Point2D v1(line4.direction[0], line4.direction[1]);
                    Point2D v2(line4.direction[2], line4.direction[3]);
                    RealType vnorm = (v2 - v1).euclidean_norm();
                    Point2D segment_coords = scaling_range * model_base_ / vnorm;

                    typename Space4D::Segment4D segment1(line4,  segment_coords);
                    typename Space4D::Segment4D segment2(line4, -segment_coords);

                    project_segment(segment1, image1, image2);
                    project_segment(segment2, image1, image2);
                }
            }
        }
    }

    void project_segment(const typename Space4D::Segment4D &segment, bo::RawImage2D<RealType> &image1,
                         bo::RawImage2D<RealType> &image2)
    {
        const RealType delta_t = RealType(0.01);

        RealType t1 = segment.second[0];
        RealType t2 = segment.second[1];

        typename Space4D::Point4D p = segment.first.point;
        typename Space4D::Point4D v = segment.first.direction;

        if (t1 > t2)
            std::swap(t1, t2);

        Vector<int, 4> tmp(0, 0, 0, 0);

        for (RealType t = t1; t < t2; t += delta_t)
        {
            typename Space4D::Point4D xt = p + v * t;

            Vector<int, 4> rounded((int)xt[0], (int)xt[1], (int)xt[2], (int)xt[3]);

            if (rounded != tmp)
            {
                tmp = rounded;

                if (rounded[0] >= 0 && rounded[0] < image1.width() &&
                    rounded[1] >= 0 && rounded[1] < image1.height())
                {
                    image1(rounded[0], rounded[1]) += 1;
                }

                if (rounded[2] >= 0 && rounded[2] < image2.width() &&
                    rounded[3] >= 0 && rounded[3] < image2.height())
                {
                    image2(rounded[2], rounded[3]) += 1;
                }
            }
        }
    }

    // Detects the references that define probable poses of the model within the 
    // given features.
    ReferenceVotes fast_detect(const Features &object_features, RealType probability, 
                               typename Space4D::Size4D divisions_per_dimension,
                               std::size_t maximal_resolution_level,
                               std::size_t cell_resolution_increment,
                               SearchArea reference_box1,
                               SearchArea reference_box2,
                               Point2D scaling_range = Point2D(0.95f, 1.05f))
    {
        ReferenceVotes ref_votes;

        // Create the root space object.
        typename Space4D::Point4D p1(reference_box1.first[0], reference_box1.first[1],
                                  reference_box2.first[0], reference_box2.first[1]);
        typename Space4D::Point4D p2(reference_box1.second[0], reference_box1.second[1],
                                  reference_box2.second[0], reference_box2.second[1]);

        Space4D s (typename Space4D::Box4D(p1, p2), divisions_per_dimension,
                   maximal_resolution_level, cell_resolution_increment, 0);

        // In the case if the scaling is incorrect.
        normalize_scaling_range(scaling_range);

        // Hierarchical search for the vote peak in the space.
        process_space(s, object_features, scaling_range, probability);

        // Get all subspaces from the last resolution level;
        typename Space4D::Spaces leafs;
        get_resolution_level(s, leafs, maximal_resolution_level);

        // Sort in descending order.
        std::sort(leafs.rbegin(), leafs.rend());

        // Extract the references.
        for (typename Space4D::Spaces::const_iterator it = leafs.begin();
             it != leafs.end(); ++it)
        {
            // Attention: the space is approximated by its mass center!
            typename Space4D::Point4D c = it->get_mass_center();

            Reference ref(Point2D(c[0], c[1]), Point2D(c[2], c[3]));
            ReferenceVote rv(ref, it->get_votes());
            ref_votes.push_back(rv);
        }

        return ref_votes;
    }

    // Detects the references that define probable poses of the model within the
    // given features.
    ReferenceVotes cross_level_detect(const Features &object_features, RealType probability,
                                      typename Space4D::Size4D divisions_per_dimension,
                                      std::size_t maximal_resolution_level,
                                      std::size_t cell_resolution_increment,
                                      SearchArea reference_box1,
                                      SearchArea reference_box2,
                                      Point2D scaling_range = Point2D(0.95f, 1.05f))
    {
        ReferenceVotes ref_votes;

        // Create the root space object.
        typename Space4D::Point4D p1(reference_box1.first[0], reference_box1.first[1],
                                  reference_box2.first[0], reference_box2.first[1]);
        typename Space4D::Point4D p2(reference_box1.second[0], reference_box1.second[1],
                                  reference_box2.second[0], reference_box2.second[1]);

        Space4D s (typename Space4D::Box4D(p1, p2), divisions_per_dimension,
                   maximal_resolution_level, cell_resolution_increment, 0);

        // In the case if the scaling is incorrect.
        normalize_scaling_range(scaling_range);

        std::size_t divisions = divisions_per_dimension[0] * divisions_per_dimension[1] *
                divisions_per_dimension[2] * divisions_per_dimension[3];

        // Initialize temporary collections of spaces.
        typename Space4D::Spaces in, out;
        in.push_back(s);

        for (std::size_t level = 0; level < maximal_resolution_level; ++level)
        {
            out.clear();
            out.reserve(in.size() * divisions);

            // Compute the votes for the level subdivision spaces.
            for (typename Space4D::Spaces::iterator it = in.begin(); it != in.end(); ++it)
            {
                subdivide_with_votes(*it, object_features, scaling_range);
                // Insert all subspaces into the output collection.
                out.insert(out.end(), it->get_subspaces().begin(), it->get_subspaces().end());
            }

            // Sort the subspaces in descending order.
            std::sort(out.rbegin(), out.rend());

            // Find minimal number of subspaces that satisfy the probability constrains.
            std::size_t n;
            if (level == maximal_resolution_level - 1)
            {
                n = out.size();
            }
            else
            {
                n = SubPolicy::probabilistic(out, probability);
            }

            // Update the input collection with the best subspaces.
            in.clear();
            in.insert(in.end(), out.begin(), out.begin() + n);
        }

        // Extract the references.
        for (typename Space4D::Spaces::const_iterator it = in.begin(); it != in.end(); ++it)
        {
            // Attention: the space is approximated by its mass center!
            typename Space4D::Point4D c = it->get_mass_center();

            Reference ref(Point2D(c[0], c[1]), Point2D(c[2], c[3]));
            ReferenceVote rv(ref, it->get_votes());
            ref_votes.push_back(rv);
        }

        return ref_votes;
    }


    // Reconstructs the points of the model that has the pose defined by the given reference.
    Points2D reconstruct(const Reference &reference_points)
    {
        Points2D points;

        // Reconstruct the model points from all alpha and beta angles and
        // the given two reference points.
        for (typename DualPointGHT<RealType>::ATable::const_iterator ait = atable_.begin();
             ait != atable_.end(); ++ait)
            for (typename DualPointGHT<RealType>::ATableRow::const_iterator it = ait->begin();
                 it != ait->end(); ++it)
            {
                points.push_back(find_intersection(reference_points, *it));
            }

        return points;
    }

private:

    Reference model_reference_;
    RealType model_base_;
    RealType tangent_accuracy_;
    ATable atable_;
    RealType pi_;

    // Row index in the alpha-table for the given tangent angle.
    inline std::size_t atable_index(RealType gamma)
    {
        return static_cast<std::size_t>((gamma + pi_) / (2 * pi_) * (atable_.size() - 1));
    }

    // Approximate gamma angle that corresponds to the alpha-table row index.
    inline RealType atable_gamma(std::size_t index)
    {
        return index * 2 * pi_ / (atable_.size() - 1) - pi_;
    }

    // Rotates the given vector in a radians.
    inline Point2D rotate(const Point2D &v, RealType a) const
    {
        RealType cosa = std::cos(a);
        RealType sina = std::sin(a);

        return Point2D(cosa * v[0] - sina * v[1], sina * v[0] + cosa * v[1]);
    }

    // Computes the "positive" normal to the given vector.
    inline Point2D normal(const Point2D &v)
    {
        return rotate(v, pi_ / 2);
    }

    // Computes the signed angle in radians [-pi, pi] betwen the vectors relatively to
    // the base vector.
    RealType angle(const Point2D &base, const Point2D &v)
    {
        // Cosine between the vectors.
        RealType cosa = base * v / base.euclidean_norm() / v.euclidean_norm();

        // Angle without the sign.
        RealType a = std::acos(cosa);

        // Normal vector for the base.
        Point2D norm_base = normal(base);

        // The sign is defined by the halfspace relatively to base where v is located.
        int sign = norm_base * v < 0 ? -1 : 1;

        // Angle with sign.
        return sign * a;
    }

    // Computes the intersection of two lines, defined by the reference points and
    // the given element of the alpha-table.
    Point2D find_intersection(const Reference &reference_points, const ATableElement &e)
    {
        // Reference points.
        Point2D p1 = reference_points.first;
        Point2D p2 = reference_points.second;

        // Reference vector.
        Point2D ab = p2 - p1;

        // If the intersection is located on the reference line, compute the intersection
        // directly.
        if (e.first == 0)
        {
            return p1 + e.second * ab;
        }

        // Find the first direction with the base norm.
        Point2D v1 = rotate(ab, e.first);

        // Calculate the intersection point.
        RealType sinb = std::sin(e.second);
        RealType sinba = std::sin(e.second - e.first);

        // The intersection is not on the reference line and not in the infinity.
        BOOST_ASSERT(sinba != 0);

        // The coordinate of the intersection point relatively to v1.
        RealType t = sinb / sinba;

        return p1 + t * v1;
    }


    // Fills in the alpha-table using the given model features
    // and two reference points.
    void encode(const Features &model_features)
    {
        const RealType epsilon = RealType(0.001);

        // Define the number of discrete tangent angles.
        unsigned int tangent_angle_number = static_cast<unsigned int>(2 * pi_ / tangent_accuracy_);

        // Allocate memory for the alpha-table.
        atable_.resize(tangent_angle_number);

        // The reference vector.
        Point2D ab = model_reference_.second - model_reference_.first;

        // Fill in the alpha-table.
        for (typename Features::const_iterator it = model_features.begin();
             it != model_features.end(); ++it)
        {
            // Current model point.
            Point2D c = it->first;
            // Boundary tangent at the current model point.
            Point2D tangent = it->second;

            // Vectors from the reference points to the current model point.
            Point2D v1 = c - model_reference_.first;
            Point2D v2 = c - model_reference_.second;

            // Compute the reference angles.
            RealType alpha = angle(ab, v1);
            RealType beta = angle(ab, v2);
            // Compute the tangential angle.
            RealType gamma = angle(v1, tangent);

            // Correction for the points located on the reference line.
            // This case is encoded as: alpha = 0; beta = coordinate of the current model point
            // on the reference line relatively to the reference vector.
            if (std::abs(std::sin(beta)) < epsilon)
            {
                alpha = 0;
                RealType abnorm =  ab.euclidean_norm();
                beta = ab * v1 / (abnorm * abnorm);
            }

            // Define the row for the computed angles in the alpha-table and insert
            // the angles into the table.
            atable_.at(atable_index(gamma)).push_back(ATableElement(alpha, beta));
        }
    }

    void normalize_scaling_range(Point2D &scaling_range)
    {
        if (scaling_range[0] < 0)
            scaling_range[0] = 0;

        if (scaling_range[1] < 0)
            scaling_range[1] = 0;

        if (scaling_range[0] > scaling_range[1])
            std::swap(scaling_range[0], scaling_range[1]);
    }

    inline void subdivide_with_votes(Space4D &s, const Features &object_features,
                                     const Point2D &scaling_range)
    {
        // Create the space subdivision.
        s.subdivide();

        // Calculate the votes for the obtained subspaces.
        for (typename Space4D::Spaces::iterator it = s.get_subspaces().begin();
             it != s.get_subspaces().end(); ++it)
        {
            feature_to_vote(*it, object_features, scaling_range);
        }
    }

    // Recursively fills in the space tree calculating votes for each subspace.
    void process_space(Space4D &s, const Features &object_features, const Point2D &scaling_range,
                       RealType probability)
    {
        if (s.get_resolution_level() >= s.get_max_resolution_level())
            return;

        // Subdivide the space and compute votes for its subspaces.
        subdivide_with_votes(s, object_features, scaling_range);

        // Sorting subspaces in descending order!
        std::sort(s.get_subspaces().rbegin(), s.get_subspaces().rend());

        // The minimal number of spaces from the beginning of the space collection such
        // that the probability of maximal element is not less than the given value.
        std::size_t n = SubPolicy::probabilistic(s.get_subspaces(), probability);

        for (std::size_t i = 0; i < n; ++i)
        {
            // Continue the subdivision procedure recursively.
            process_space(s.get_subspaces().at(i), object_features, scaling_range, probability);
        }
    }

    // Intersects the given space with the lines produced by the object features and increase the
    // number of the space votes.
    void feature_to_vote(Space4D &s, const Features &object_features, const Point2D &scaling_range)
    {
        for (typename Features::const_iterator it = object_features.begin(); it != object_features.end(); ++it)
        {
            feature_to_vote2(s, *it, scaling_range);
        }
    }

    // Intersects the given space with the line defined by the feature (position and tangent)
    // and increase the number of the space votes.
    inline void feature_to_vote(Space4D &s, const Feature &f, const Point2D &scaling_range)
    {
        // Reconstruct all the 4D lines from the alpha-table relatively to the current feature.
        for (std::size_t index = 0; index < atable_.size(); ++index)
        {
            // Probable angle between the directional vector and the tangent.
            RealType gamma = atable_gamma(index);

            // Find all corresponding directional vectors v2 defined by the (alpha, beta) angles of
            // the current table row.
            for (typename ATableRow::const_iterator abit = atable_.at(index).begin();
                 abit !=  atable_.at(index).end(); ++abit)
            {
                typename Space4D::Line4D line4 = line4_from_feature_and_atable_element(f, gamma, *abit);

                // Create two segments on this line that correspond to the given scaling range.
                Point2D v1(line4.direction[0], line4.direction[1]);
                Point2D v2(line4.direction[2], line4.direction[3]);
                RealType vnorm = (v2 - v1).euclidean_norm();
                Point2D segment_coords = scaling_range * model_base_ / vnorm;

                typename Space4D::Segment4D segment1(line4,  segment_coords);
                typename Space4D::Segment4D segment2(line4, -segment_coords);

                // Intersect and compute the votes.
                s.vote_maxnorm(s.intersect(segment1));
                s.vote_maxnorm(s.intersect(segment2));
            }
        }
    }

    inline typename Space4D::Line4D line4_from_feature_and_atable_element(const Feature &f,
                                                                          RealType gamma,
                                                                          const ATableElement &e)
    {
        Point2D c = f.first;
        Point2D tan = f.second;

        // Normalized directional vector.
        Point2D v1 = rotate(tan, -gamma);
        v1 = v1 / v1.euclidean_norm();

        RealType alpha = e.first;
        RealType beta = e.second;

        Point2D v2 = rotate(v1, beta - alpha);

        // Consider the directional vectors ratio.
        if (alpha != 0)
            // The conventional case.
            v2 *= std::sin(alpha) / std::sin(beta);
        else
            // The directional vectors are located on the reference line.
            v2 *= 1 - 1 / beta;

        // Compose a 4D line.
        // The point on this line.
        typename Space4D::Point4D p4(c.x(), c.y(), c.x(), c.y());
        // The directional vector of the line.
        typename Space4D::Point4D v4(v1.x(), v1.y(), v2.x(), v2.y());
        // The line in 4D.
        typename Space4D::Line4D line4(p4, v4);

        return line4;
    }

    // Intersects the given space with the hyperplane defined by the feature (position and tangent)
    // and increase the number of the space votes.
    inline void feature_to_vote2(Space4D &s, const Feature &f, const Point2D &scaling_range)
    {
        // Reconstruct all the 4D hyperplanes from the alpha-table relatively to the current feature.
        for (std::size_t index = 0; index < atable_.size(); ++index)
        {
            // Probable angle between the directional vector and the tangent.
            RealType gamma = atable_gamma(index);

            // Find all corresponding directional vectors v2 defined by the (alpha, beta) angles of
            // the current table row.
            for (typename ATableRow::const_iterator abit = atable_.at(index).begin();
                 abit !=  atable_.at(index).end(); ++abit)
            {
                // Reconstruct the 4D line.
                typename Space4D::Line4D line4 = line4_from_feature_and_atable_element(f, gamma, *abit);

                // Use the hyperplane paradigm.
                typename Space4D::Point4D norm(line4.direction[3], line4.direction[2],
                                            -line4.direction[1], -line4.direction[0]);
                typename Space4D::Point4D point = line4.point;
                typename Space4D::Hyperplane plane(point, norm);

                // Check the plane constrains here and vote.
                if (is_in_scaling_constrains(s, plane, scaling_range) &&
                    is_in_tangent_constrains(s, plane, RealType(0.3)))
                {
                     s.vote_descrete(plane);
                }
            }
        }
    }

    int location_by_line(const Point2D &bot, const Point2D &top,
                         const Point2D &b, const Point2D &normal) const
    {
        Point2D d = top - bot;

        RealType proj_min = std::numeric_limits<RealType>::max();
        RealType proj_max = -proj_min;

        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
            {
                Point2D p = bot + Point2D(d[0] * i, d[1] * j);

                RealType proj = (p - b) * normal;

                if (proj_min > proj)
                    proj_min = proj;

                if (proj_max < proj)
                    proj_max = proj;
            }

        if ((proj_min == 0 && proj_max == 0) || proj_min * proj_max < 0)
            return 0;

        return proj_max > 0 ? 1 : -1;
    }

    bool is_in_scaling_constrains(const Space4D &s,
                                  const typename Space4D::Hyperplane &plane,
                                  const Point2D &scaling_range) const
    {
        // Compute directional vectors.
        Point2D v1(-plane.normal()[3], -plane.normal()[2]);
        Point2D v2(plane.normal()[1], plane.normal()[0]);

        RealType vnorm = (v2 - v1).euclidean_norm();
        Point2D scaling = scaling_range * model_base_ / vnorm;

        typename Space4D::Box4D box = s.get_bounding_box();

        // Box projections.
        Point2D bot_v1(box.first[0], box.first[1]);
        Point2D top_v1(box.second[0], box.second[1]);
        Point2D bot_v2(box.first[2], box.first[3]);
        Point2D top_v2(box.second[2], box.second[3]);

        Point2D b(plane.point()[0], plane.point()[1]);

        // Points on the lines.
        Point2D b_v1q1 = b + v1 * scaling[0];
        Point2D b_v1q2 = b + v1 * scaling[1];
        Point2D b_v2q1 = b + v2 * scaling[0];
        Point2D b_v2q2 = b + v2 * scaling[1];

        Point2D b_v1q1_inv = b - v1 * scaling[0];
        Point2D b_v1q2_inv = b - v1 * scaling[1];
        Point2D b_v2q1_inv = b - v2 * scaling[0];
        Point2D b_v2q2_inv = b - v2 * scaling[1];

        return (location_by_line(bot_v1, top_v1, b_v1q1, v1) >= 0 &&
                location_by_line(bot_v1, top_v1, b_v1q2, v1) <= 0 &&
                location_by_line(bot_v2, top_v2, b_v2q1, v2) >= 0 &&
                location_by_line(bot_v2, top_v2, b_v2q2, v2) <= 0) ||
               (location_by_line(bot_v1, top_v1, b_v1q1_inv, -v1) >= 0 &&
                location_by_line(bot_v1, top_v1, b_v1q2_inv, -v1) <= 0 &&
                location_by_line(bot_v2, top_v2, b_v2q1_inv, -v2) >= 0 &&
                location_by_line(bot_v2, top_v2, b_v2q2_inv, -v2) <= 0);
    }

    inline void sector_constraint(const Point2D &v, RealType sigma,
                                  Point2D &nw) const
    {
        Point2D w = rotate(v, sigma);

        nw[0] = -w[1];
        nw[1] = w[0];

        // Correct the sign of the normal.
        if (nw * (v - w) > 0)
            nw *= -1;
    }

    bool is_in_tangent_constrains(const Space4D &s,
                                  const typename Space4D::Hyperplane &plane,
                                  const RealType sigma) const
    {
        // Compute directional vectors v1 and v2.
        Point2D v1(-plane.normal()[3], -plane.normal()[2]);
        Point2D v2(plane.normal()[1], plane.normal()[0]);

        typename Space4D::Box4D box = s.get_bounding_box();

        // Box projections.
        Point2D bot_v1(box.first[0], box.first[1]);
        Point2D top_v1(box.second[0], box.second[1]);
        Point2D bot_v2(box.first[2], box.first[3]);
        Point2D top_v2(box.second[2], box.second[3]);

        Point2D b(plane.point()[0], plane.point()[1]);

        // Calculate vectors constraining the sigma-sector.
        Point2D nw1, nw2, nu1, nu2;
        sector_constraint(v1, sigma, nw1);
        sector_constraint(v1, -sigma, nu1);
        sector_constraint(v2, sigma, nw2);
        sector_constraint(v2, -sigma, nu2);

        return (location_by_line(bot_v1, top_v1, b, nw1) *
                location_by_line(bot_v1, top_v1, b, nu1) >= 0) &&
               (location_by_line(bot_v2, top_v2, b, nw2) *
                location_by_line(bot_v2, top_v2, b, nu2) >= 0);
    }

    // Recursively traces the space tree and inserts into the container the elements from
    // the given resolution level.
    void get_resolution_level(const Space4D &s, typename Space4D::Spaces &container,
                              std::size_t resolution_level)
    {
        if (s.get_resolution_level() == resolution_level)
            container.push_back(s);
        else
        {
            const typename Space4D::Spaces subs = s.get_subspaces();

            for (typename Space4D::Spaces::const_iterator it = subs.begin(); it != subs.end(); ++it)
            {
                get_resolution_level(*it, container, resolution_level);
            }
        }

    }

};

} // namespace recognition
} // namespace methods
} // namespace bo

#endif // DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419_
