
/******************************************************************************

  parallel_planes_tiling.hpp, v 1.0.0 2012.11.27

  Implementation of several surface tiling methods, working with parallel planes.

  Copyright (c) 2012
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

#ifndef PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
#define PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "bo/config.hpp"
#include "bo/mesh.hpp"

namespace bo {
namespace methods {
namespace surfaces {

template <typename RealType>
class MaxSpanPropagation: public boost::noncopyable
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> ParallelPlane;
    typedef Mesh<RealType> Mesh;
    typedef boost::shared_ptr<ParallelPlane> ParallelPlanePtr;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
