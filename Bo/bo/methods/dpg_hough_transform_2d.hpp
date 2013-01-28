
/******************************************************************************

  dpg_hough_transform_2d.hpp, v 1.0.1 2013.01.12

  Fast dual-point generalized Hough Transform for 2D object recognition.

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
#include <cassert>

#include <boost/noncopyable.hpp>
#include <boost/math/constants/constants.hpp>

#include "bo/vector.hpp"
#include "bo/blas/blas.hpp"

#include <iostream>

namespace bo {
namespace methods {
namespace recognition {

namespace detail{

template <typename RealType>
class Space
{
public:

    typedef Vector<RealType, 4> Point4D;
    typedef Vector<int, 4> Index4D;
    typedef std::pair<Point4D, Point4D> Box4D;
    // A line in 4D is defined by any point located
    // on this line and its directional vector.
    typedef std::pair<Point4D, Point4D> Line4D;
    typedef Vector<RealType, 2> SegmentCoordinates;
    typedef std::vector<Space> Spaces;
    // A segment in 4D is modelled by a 4D line and the
    // coordinates of two points on this line
    // relatively to the line's direction vector.
    typedef std::pair<Line4D, SegmentCoordinates> Segment4D;

    Space(const Box4D &box = Box4D(Point4D(0,0,0,0), Point4D(0,0,0,0)), std::size_t resolution_level = 0):
        box_(box), resolution_level_(resolution_level), votes_(0)
    {
    }

    // Subdivides the space.
    void subdivide(std::size_t divisions_per_dimension)
    {
        subspaces_.reserve(divisions_per_dimension * divisions_per_dimension *
                           divisions_per_dimension * divisions_per_dimension);

        Point4D sizes4d = box_.second - box_.first;
        Point4D steps4d = sizes4d / divisions_per_dimension;

        for (std::size_t d0 = 0; d0 < divisions_per_dimension; ++d0)
            for (std::size_t d1 = 0; d1 < divisions_per_dimension; ++d1)
                for (std::size_t d2 = 0; d2 < divisions_per_dimension; ++d2)
                    for (std::size_t d3 = 0; d3 < divisions_per_dimension; ++d3)
                    {
                        // Create a subspace.
                        Point4D translation(d0 * steps4d[0], d1 * steps4d[1], d2 * steps4d[2], d3 * steps4d[3]);
                        Box4D b(box_.first + translation, box_.first + translation + steps4d);
                        Space s(b, resolution_level_ + 1);

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

    void reset_votes()
    {
        votes_ = 0;
    }

    inline std::size_t get_resolution_level() const
    {
        return resolution_level_;
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
        Point4D s = segment.first.second * (segment.second[1] - segment.second[0]);

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

    // Increase the number of votes on the number of intersection of the given segment
    // with the grid of the given cell size in the current space.
    // Attention: the space must contain integer number of cells (in each dimension)!
    void vote_descrete(const Segment4D &segment, const Point4D &cell_size)
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
        Point4D v = segment.first.second;

        // For each dimension define the coordinate step beetween two cells (dimension levels).
        for (std::size_t d = 0; d < 4; ++d)
        {
            if (v[d] != 0)
            // If it is zero, the segment does not intersect the levels of this dimension and
            // we can skip further analysis.
            {
                RealType dt = std::abs(cell_size[d] / v[d]);

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

    // The number of votes is equal to the taxicab norm of the corresponding inner
    // segment relatively to the grid with the given cell size.
    void vote_taxicab(const Segment4D &segment, const Point4D &cell_size)
    {
        // If the segment is not zero.
        if (segment.second[0] == segment.second[1])
            return;

        std::size_t taxicab_dst = 0;

        // Local coordinate origin.
        Point4D bot = box_.first;
        Point4D top = box_.second;

        // The end points of the segment.
        Point4D p1 = segment.first.first + segment.first.second * segment.second[0];
        Point4D p2 = segment.first.first + segment.first.second * segment.second[1];

        // Segment direction.
        Point4D v = segment.first.second;

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

                assert (proj1 >= bot[d] && proj2 >= bot[d]);

                RealType block1 = std::floor((proj1 - bot[d]) / cell_size[d]);
                RealType block2 =  std::ceil((proj2 - bot[d]) / cell_size[d]);

                taxicab_dst += static_cast<std::size_t>(block2 - block1);
            }
        }

        votes_ += taxicab_dst;
    }

    Segment4D intersect(Segment4D segment)
    {
        // Cut the segment in each of four dimensions.
        for (std::size_t d = 0; d < 4; ++d)
            cut_segment(segment, d, box_.first[d], box_.second[d]);

        return segment;
    }


    Segment4D intersect(Line4D line)
    {
        // Create the "infinite" segment coordinates that models the whole line.
        RealType mmax = std::numeric_limits<RealType>::max();
        SegmentCoordinates coord(-mmax, mmax);

        // Create the "infinite" segment.
        Segment4D infinite_seg(line, coord);

        return intersect(infinite_seg);
    }


    // Operator is needed for sorting.
    bool operator < (const Space &other) const
    {
        return votes_ < other.get_votes() ? true : false;
    }

private:

    Spaces subspaces_;
    Box4D box_;
    std::size_t resolution_level_;
    RealType votes_;


    // Cuts the segment in the given dimension according to two given levels.
    void cut_segment(Segment4D &seg, std::size_t dimension, RealType level1, RealType level2)
    {
        // Sort the levels: level1 <= level2.
        order(level1, level2);

        // The origin point and the direction.
        Point4D p = seg.first.first;
        Point4D v = seg.first.second;

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

    DualPointGHT(const Features &model_features, const Reference &model_reference, 
                 RealType tangent_accuracy = RealType(0.005)):
    model_reference_(model_reference), tangent_accuracy_(tangent_accuracy),
    pi_(boost::math::constants::pi<RealType>())
    {
        encode(model_features);

        model_base_ = (model_reference_.second - model_reference_.first).euclidean_norm();
    }

    // Detects the references that define probable poses of the model within the 
    // given features.
    ReferenceVotes fast_detect(const Features &object_features, RealType probability, 
                               unsigned int divisions_per_dimension,
                               unsigned int maximal_resolution_level,
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
        Space4D s(typename Space4D::Box4D(p1, p2));

        // Initialize the grid size.
        std::size_t divisions = std::pow(divisions_per_dimension, maximal_resolution_level + 1);
        grid_size_ = (p2 - p1) / divisions;

        // In the case if the scaling is incorrect.
        normalize_scaling_range(scaling_range);

        // Hierarchical search for the vote peak in the space.
        process_space(s, object_features, probability, divisions_per_dimension,
                      maximal_resolution_level, scaling_range);

        // Get all subspaces from the last resolution level;
        typename Space4D::Spaces leafs;
        get_resolution_level(s, leafs, maximal_resolution_level);

        std::sort(leafs.begin(), leafs.end());

        // Extract the references.
        for (typename Space4D::Spaces::const_reverse_iterator it = leafs.rbegin();
             it != leafs.rend(); ++it)
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
    typename Space4D::Point4D grid_size_;

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
    inline Point2D rotate(const Point2D &v, RealType a)
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

    // Computes the signed angle in radians [-pi/2, pi/2] betwen the vectors relatively to
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

        // Find the first normalized direction.
        Point2D v1 = rotate(ab, e.first);
        v1 = v1 / v1.euclidean_norm();

        // Find the second normalized direction.
        Point2D v2 = rotate(ab, e.second);
        v2 = v2 / v2.euclidean_norm();

        // Calculate the intersection point.
        RealType sina = std::sin(e.first);
        RealType sinb = std::sin(e.second);

        // The intersection is not on the reference line.
        BOOST_ASSERT(sinb != 0);

        // The coordinate of the intersection point relatively to v1.
        RealType t = ab[0] / (v1[0] - v2[0] * sina / sinb);

        return p1 + t * v1;
    }


    // Fills in the alpha-table using the given model features
    // and two reference points.
    void encode(const Features &model_features)
    {
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
            if (std::sin(beta) == 0)
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

    // Recursively fills in the space tree calculating votes for each subspace.
    void process_space(Space4D &s, const Features &object_features, RealType probability,
                       RealType divisions_per_dimension, RealType maximal_resolution_level,
                       const Point2D &scaling_range)
    {
        if (s.get_resolution_level() >= maximal_resolution_level)
            return;

        // Create the space subdivision.
        s.subdivide(divisions_per_dimension);
        //typename Space4D::Spaces subs = s.get_subspaces();

        // Calculate the votes for the obtained subspaces.
        for (typename Space4D::Spaces::iterator it = s.get_subspaces().begin();
             it != s.get_subspaces().end(); ++it)
        {
            feature_to_vote(*it, object_features, scaling_range);
        }

        // Subspaces with less votes first.
        std::sort(s.get_subspaces().begin(), s.get_subspaces().end());

        // Continue the subdivision procedure for the ONE subspace with the maximal number of votes.
        process_space(s.get_subspaces().back(), object_features, probability, divisions_per_dimension,
                      maximal_resolution_level, scaling_range);
        /*
        for (typename Space4D::Spaces::reverse_iterator it = s.get_subspaces().rbegin(); it != s.get_subspaces().rend(); ++it)
        {
            process_space(*it, object_features, probability, divisions_per_dimension,
                          maximal_resolution_level);
        }
        */
    }

    // Intersects the given space with the lines produced by the object features and increase the
    // number of the space votes.
    void feature_to_vote(Space4D &s, const Features &object_features, const Point2D &scaling_range)
    {
        for (typename Features::const_iterator it = object_features.begin(); it != object_features.end(); ++it)
        {
            feature_to_vote(s, *it, scaling_range);
        }
    }

    // Intersects the given space with the line defined by the feature (position and tangent)
    // and increase the number of the space votes.
    inline void feature_to_vote(Space4D &s, const Feature &f, const Point2D &scaling_range)
    {
        Point2D c = f.first;
        Point2D tan = f.second;

        // Reconstruct all the 4D lines from the alpha-table relatively to the current feature.
        for (std::size_t index = 0; index < atable_.size(); ++index)
        {
            // Probable angle between the directional vector and the tangent.
            RealType gamma = atable_gamma(index);

            // Probable normalized directional vector.
            Point2D v1 = rotate(tan, -gamma);
            v1 = v1 / v1.euclidean_norm();

            // Find all corresponding directional vectors v2 defined by the (alpha, beta) angles of
            // the current table row.
            for (typename ATableRow::const_iterator abit = atable_.at(index).begin();
                 abit !=  atable_.at(index).end(); ++abit)
            {
                RealType alpha = abit->first;
                RealType beta = abit->second;

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

                // Create two segments on this line that correspond to the given scaling range.
                RealType vnorm = (v2 - v1).euclidean_norm();
                Point2D segment_coords = scaling_range * model_base_ / vnorm;

                typename Space4D::Segment4D segment1(line4,  segment_coords);
                typename Space4D::Segment4D segment2(line4, -segment_coords);

                // Intersect and compute the votes.
                s.vote_taxicab(s.intersect(segment1), grid_size_);
                s.vote_taxicab(s.intersect(segment2), grid_size_);
            }
        }
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
