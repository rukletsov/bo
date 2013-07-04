
/******************************************************************************

  Helper types used in surface reconstruction routines.

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

#ifndef TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021
#define TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021

#include <vector>

#include "bo/core/vector.hpp"
#include "bo/math/pca.hpp"
#include "bo/math/mean.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// Represents a point cloud located inside a thin disk. Though the disk can be
// represented by a median plane, slight deviations from it may exist.
template <typename RealType>
class PointsDisk3D
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> PlaneData;

    PointsDisk3D(const PlaneData& data): data_(data)
    {
        // Compute plane origin.
        origin_ = bo::math::mean(data_);

        // Employ PCA to estimate plane normal.
        typedef math::PCA<RealType, 3> PCAEngine;
        PCAEngine pca;
        typename PCAEngine::Result result = pca(data_);
        normal_ = result.template get<1>()[0];
    }

    const PlaneData& data() const
    {
        return data_;
    }

    Point3D origin() const
    {
        return origin_;
    }

    Point3D normal() const
    {
        return normal_;
    }

private:
    PlaneData data_;
    Point3D origin_;
    Point3D normal_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021
