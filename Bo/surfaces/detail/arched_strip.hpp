
/******************************************************************************

  Functor that detects whether a point lies inside an arched strip of a
  certain form with additional conditions.

  Copyright (c) 2013
  Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef ARCHED_STRIP_HPP_729C21F2_5247_4198_A119_257F925B6BD2
#define ARCHED_STRIP_HPP_729C21F2_5247_4198_A119_257F925B6BD2

#include <boost/function.hpp>

#include "bo/core/vector.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// Represents an arched strip with additional conditions. Call operator() to determine
// whether a point lies inside the strip.
template <typename RealType, std::size_t Dim>
class ArchedStrip
{
public:
    typedef Vector<RealType, Dim> Point;
    typedef boost::function<RealType (Point, Point)> Metric;

public:
    ArchedStrip(Point ref_pt): ref_pt_(ref_pt)
    { }

    ArchedStrip(Point ref_pt, RealType delta_min, RealType delta_max, Point prop,
                Metric dist_fun):
        ref_pt_(ref_pt), delta_min_(delta_min), delta_max_(delta_max), prop_(prop),
        dist_fun_(dist_fun)
    { }

    bool operator()(const Point& pt)
    {
        // We need to check three conditions are met:
        // 1. the point lies outside the delta_min circle;
        // 2. the point lies inside the delta_max circle;
        // 3. the point lies in front of the diving plane normal to
        //    propagation vector.

        RealType dist = dist_fun_(ref_pt_, pt);
        bool cond1 = (dist >= delta_min_);
        bool cond2 = (dist <= delta_max_);

        // Dot product is positive if the angle between vectors < 90 deg.
        bool cond3 = ((prop_ * (pt - ref_pt_)) > 0);

        return (cond1 && cond2 && cond3);
    }

private:
    Point ref_pt_;
    RealType delta_min_;
    RealType delta_max_;
    Point prop_;
    Metric dist_fun_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // ARCHED_STRIP_HPP_729C21F2_5247_4198_A119_257F925B6BD2
