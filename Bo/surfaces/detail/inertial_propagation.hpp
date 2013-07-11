
/******************************************************************************

  Helper types related to the inertial component of the complex propagation.

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

#ifndef INERTIAL_PROPAGATION_HPP_FEAC7E0E_EA26_11E2_94AD_DD539A890B3C
#define INERTIAL_PROPAGATION_HPP_FEAC7E0E_EA26_11E2_94AD_DD539A890B3C

#include <boost/shared_ptr.hpp>

#include "bo/core/vector.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// Standard inertial propagation that uses provided point as a previous one.
template <typename RealType>
class InertialPropagation
{
public:
    typedef InertialPropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    InertialPropagation(RealType weight): weight_(weight)
    { }

    virtual Point3D get(const Point3D& current, const Point3D& previous) const
    {
        Point3D inertial = current - previous;

        // TODO: remove this by redesigning the algorithm and requiring inertial
        // vector to be non-zero.
        Point3D inertial_normalized(0);
        try
        {
            // TODO: remove this block by refactoring bo::Vector class.
            // Normalize vector.
            inertial_normalized.assign(inertial.normalized(), 3);
        }
        catch (...)
        { }

        return weight_ * inertial_normalized;
    }

protected:
    RealType weight_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // INERTIAL_PROPAGATION_HPP_FEAC7E0E_EA26_11E2_94AD_DD539A890B3C
