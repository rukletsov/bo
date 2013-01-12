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

#include "bo/vector.hpp"

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
    typedef std::vector<std::list<ATableElement> > ATable;
     

    DualPointGHT(const Features &model_features, const Reference &model_reference):
                model_reference_(model_reference)
    {
        encode(model_features, model_reference);
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

        return points;
    }

private:

    Reference model_reference_;
    ATable atable_;

    // Fills in the alpha-table using the given model features
    // and two reference points.
    void encode(const Features &model_features, const Reference &model_reference)
    {
        
    }

};

} // namespace recognition.

} // namespace methods.
} // namespace bo.

#endif //DPG_HOUGH_TRANSFORM_HPP_55B0255B_E6C8_4302_9114_D1B684CD3419
