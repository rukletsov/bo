
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

#ifndef DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419
#define DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419

#include <utility>
#include <vector>
#include <list>
#include <cmath>
#include <boost/noncopyable.hpp>

#include "bo/vector.hpp"
#include "bo/blas/blas.hpp"

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
    model_reference_(model_reference), tangent_accuracy_(tangent_accuracy)
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
        for (ATable::const_iterator ait = atable_.begin(); ait != atable_.end(); ++ait)
            for (ATableRow::const_iterator it = ait->begin(); it != ait->end(); ++it)
            {
                points.push_back(find_intersection(reference_points, *it)); 
            }

        return points;
    }

private:

    Reference model_reference_;
    RealType tangent_accuracy_;
    ATable atable_;

    // Row index in the alpha-table for the given cosine of the tangent angle.
    std::size_t atable_index(RealType cos_gamma)
    {
        return static_cast<std::size_t>((cos_gamma + 1) / 2 * (atable_.size() - 1));
    }

    // Fill in the given 2x2 matrix as a rotational one accordingly to the
    // given rotation angle cosine.
    void rotational(blas::matrix<RealType> &m, RealType cosa)
    {
        RealType sina = std::sqrt(1 - cosa * cosa);
        m(0, 0) = cosa;
        m(0, 1) = - sina; //! or minus.
        m(1, 0) = sina;
        m(1, 1) = cosa;
    }

    // Computes the intersection of two lines, defined by the reference points and 
    // the given element of the alpha-table.
    Point2D find_intersection(const Reference &reference_points, const ATableElement &e)
    {
        Point2D p;

        // Reference points.
        Point2D p1 = reference_points.first;
        Point2D p2 = reference_points.second;

        // Create the reference matrix.
        blas::matrix<RealType> ab(2, 1);
        ab(0, 0) = p1[1] - p2[1];
        ab(1, 0) = p2[0] - p1[0];

        // Initialize a rotational matrix.
        blas::matrix<RealType> r(2, 2);

        blas::matrix<RealType> a1 = blas::prod(rotational(r, e.first), ab);
        blas::matrix<RealType> a2 = blas::prod(rotational(r, e.second), ab);

        // Create the basis matrix.
        blas::matrix<RealType> a(2, 2);
        a(0, 0) = a1(0, 0);
        a(0, 1) = a1(1, 0);
        a(1, 0) = a2(0, 0);
        a(1, 1) = a2(1, 0);

        // Init an inverse matrix.
        blas::matrix<RealType> ia(2, 2);  
        bool is_inversible = blas::invert_matrix(a, ia);

        BOOST_ASSERT(is_inversible);

        // Compute the intersection.
        blas::matrix<RealType> b(2, 1); 
        b(0, 0) = - (a1(0, 0) * p1[0] + a1(1, 0) * p1[1]);
        b(1, 0) = - (a2(0, 0) * p2[0] + a2(1, 0) * p2[1]);

        blas::matrix<RealType> c = blas::prod(ia, b);

        p[0] = c(0, 0);
        p[1] = c(1, 0);

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

        // Normalized reference vector.
        Point2D ab = model_reference_.second - model_reference_.first;
        ab = ab.normalized();
        
        // Fill in the alpha-table.
        for (Features::const_iterator it = model_features.begin(); it != model_features.end(); ++it)
        {
            // Current model point.
            Point2D c = it->first;
            // Boundary tangent at the current model point.
            Point2D tangent = it->second;

            // Vectors from the reference points to the current model point.
            Point2D v1 = c - model_reference_.first;
            Point2D v2 = c - model_reference_.second;

            // Compute cosines of the reference angles.
            RealType cos_alpha = ab * v1 / v1.euclidean_norm();
            RealType cos_beta = ab * v2 / v2.euclidean_norm();

            // Compute cosine of the tangential angle.
            RealType cos_gamma = - v1 * tangent / v1.euclidean_norm() / tangent.euclidean_norm();

            // Define the row for the computed angles in the alpha-table and insert
            // the angles into the table.
            atable_.at(atable_index(cos_gamma)).push_back(ATableElement(cos_alpha, cos_beta));
        }
    }

};

} // namespace recognition
} // namespace methods
} // namespace bo

#endif // DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419
