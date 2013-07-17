
/******************************************************************************

  Helper intermediate type used for storing propagated contours.

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

#ifndef PROPAGATION_RESULT_HPP_44E2EB0A_EA2D_11E2_9542_AB62D8EF6274
#define PROPAGATION_RESULT_HPP_44E2EB0A_EA2D_11E2_9542_AB62D8EF6274

#include <vector>

#include "bo/core/vector.hpp"

namespace bo {
namespace surfaces {
namespace detail {

template <typename RealType>
struct PropagationResult
{
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> PropagatedContour;

    PropagationResult(): Aborted(false), HasHole(false)
    { }

    PropagationResult(bool maxsize_reached, bool hole_encountered, PropagatedContour pts):
        Aborted(maxsize_reached), HasHole(hole_encountered), Points(pts)
    { }

    bool Aborted;
    bool HasHole;
    PropagatedContour Points;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // PROPAGATION_RESULT_HPP_44E2EB0A_EA2D_11E2_9542_AB62D8EF6274

