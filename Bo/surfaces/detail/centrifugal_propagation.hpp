
/******************************************************************************

  Helper types related to the centrifugal component of the complex propagation.

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

#ifndef CENTRIFUGAL_PROPAGATION_HPP_63941A92_EA23_11E2_BC2A_2B71033B4CBB
#define CENTRIFUGAL_PROPAGATION_HPP_63941A92_EA23_11E2_BC2A_2B71033B4CBB

#include <boost/shared_ptr.hpp>

#include "bo/config.hpp"
#include "bo/core/vector.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// Base class for any centrifual propagation.
template <typename RealType>
class BaseCentrifugalPropagation
{
public:
    typedef BaseCentrifugalPropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    virtual Point3D get(const Point3D& current) const = 0;
};

// Centrifugal propagation that does nothing. Use this class if you want to disable
// centrifugal component of the complex propagation where an instance of centrifugal
// component is required.
template <typename RealType>
class EmptyCentrifugalPropagation: public BaseCentrifugalPropagation<RealType>
{
public:
    typedef Vector<RealType, 3> Point3D;

    EmptyCentrifugalPropagation(): zero_vector_(0)
    { }

    virtual Point3D get(const Point3D& current) const
    {
        BO_UNUSED(current);
        return zero_vector_;
    }

protected:
    Point3D zero_vector_;
};

// Standard centrifugal propagation that uses center of mass as a reference.
template <typename RealType>
class CentrifugalPropagation: public BaseCentrifugalPropagation<RealType>
{
public:
    typedef Vector<RealType, 3> Point3D;

    CentrifugalPropagation(const Point3D& center_of_mass, RealType weight):
        center_of_mass_(center_of_mass), weight_(weight)
    { }

    virtual Point3D get(const Point3D& current) const
    {
        // Save some processor tacts and reuse the input variable.
        Point3D inertial = current - center_of_mass_;

        // TODO: remove this by redesigning the algorithm and requiring inertial
        // vector to be non-zero.
        Point3D centrifugal_normalized(0);
        try
        {
            // TODO: remove this block by refactoring bo::Vector class.
            // Normalize vector.
            centrifugal_normalized.assign(inertial.normalized(), 3);
        }
        catch (...)
        { }

        return weight_ * centrifugal_normalized;
    }

protected:
    const Point3D center_of_mass_;
    const RealType weight_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // CENTRIFUGAL_PROPAGATION_HPP_63941A92_EA23_11E2_BC2A_2B71033B4CBB
