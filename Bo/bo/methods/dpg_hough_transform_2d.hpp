
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
#include <cmath>

#include <boost/noncopyable.hpp>
#include <boost/math/constants/constants.hpp>

#include "bo/vector.hpp"
#include "bo/blas/blas.hpp"

#include <iostream>

namespace bo {
namespace methods {
namespace recognition {

template <typename RealType>
class DualPointGHT: public boost::noncopyable
{
public:
    typedef DualPointGHT<RealType> this_type;
    typedef Vector<RealType, 2> Point2D;
    typedef std::vector<Point2D> Points2D;
    typedef std::pair<Point2D, Point2D> Reference;
    typedef std::pair<Point2D, Point2D> RecognitionArea;
    typedef std::pair<Reference, unsigned int> ReferenceVote;
    typedef std::vector<ReferenceVote> ReferenceVotes;
    typedef std::pair<Point2D, Point2D> Feature;
    typedef std::vector<Feature> Features; 
    typedef std::pair<RealType, RealType> ATableElement;
    typedef std::list<ATableElement> ATableRow;
    typedef std::vector<ATableRow> ATable;

    DualPointGHT(const Features &model_features, const Reference &model_reference, 
                 RealType tangent_accuracy = 0.5): 
    model_reference_(model_reference), tangent_accuracy_(tangent_accuracy),
    pi_(boost::math::constants::pi<RealType>())
    {
        encode(model_features);
    }

    // Detects the references that define probable poses of the model within the 
    // given features.
    ReferenceVotes fast_detect(const Features &object_features, RealType probability, 
                               unsigned int T, unsigned int M, RecognitionArea bounding_box) 
    {
        ReferenceVotes ref_votes;

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
                points.push_back(find_intersection2(reference_points, *it));
            }

        return points;
    }

private:

    Reference model_reference_;
    RealType tangent_accuracy_;
    ATable atable_;

    RealType pi_;

    // Row index in the alpha-table for the given tangent angle.
    std::size_t atable_index(RealType gamma)
    {
        return static_cast<std::size_t>((gamma + pi_) / (2 * pi_) * (atable_.size() - 1));
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
    Point2D find_intersection2(const Reference &reference_points, const ATableElement &e)
    {
        Point2D p;

        // Reference points.
        Point2D p1 = reference_points.first;
        Point2D p2 = reference_points.second;

        // Reference vector.
        Point2D ab = p2 - p1;

        // If the intersection is located on the reference line compute the intersection
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

        return p;
    }


    // Fills in the alpha-table using the given model features
    // and two reference points.
    void encode(const Features &model_features)
    {
        // Define the number of discrete tangent angles.
        unsigned int tangent_angle_number = static_cast<unsigned int>(180 / tangent_accuracy_);

        // Allocate memory for the alpha-table.
        atable_.resize(tangent_angle_number);

        // The reference vector.
        Point2D ab = model_reference_.second - model_reference_.first;

        // Fill in the alpha-table.
        for (typename DualPointGHT<RealType>::Features::const_iterator it = model_features.begin();
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

};

} // namespace recognition
} // namespace methods
} // namespace bo

#endif // DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419_
