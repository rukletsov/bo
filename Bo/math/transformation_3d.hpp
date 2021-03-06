
/******************************************************************************

  3D space transformations.

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

#ifndef TRANSFORMATION_3D_HPP_23877F3F_5EB9_4221_A422_C239EC9AA5D6_
#define TRANSFORMATION_3D_HPP_23877F3F_5EB9_4221_A422_C239EC9AA5D6_

#include <iostream>

#include "bo/core/vector.hpp"
#include "bo/math/blas_extensions.hpp"
#include "bo/math/blas_conversions.hpp"

namespace bo {
namespace math {

// A class representing transformations in 3D. Internal representation is a row-major
// 4x4 transformation matrix, but pairs (quaternion, translation) are also supported.
template <typename RealType>
class Transformation3D
{
public:
    typedef Transformation3D<RealType> self_type;
    typedef Vector<RealType, 3> Point3D;
    typedef Vector<RealType, 3> Translation;
    typedef Vector<RealType, 4> Quaternion;
    typedef math::bounded_matrix<RealType, 3, 3> RotationMatrix;
    typedef math::bounded_matrix<RealType, 4, 4> TransformationMatrix;

public:
    // Creates either identity transformation or the transformation based on the
    // given combinations of quaternion, rotation matrix, translation vector.
    Transformation3D();
    Transformation3D(const Quaternion& q);
    Transformation3D(const RotationMatrix& rot_matrix);
    Transformation3D(const Translation& t);
    Transformation3D(const Quaternion& q, const Translation& t);
    Transformation3D(const RotationMatrix& rot_matrix, const Translation& t);
    Transformation3D(const TransformationMatrix& trans_matrix);

    // Copy c-tor, d-tor and assignment operator are fine.

    // Multiplies the given 3D point on the left by the current transformation.
    Point3D operator*(Point3D point) const;

    // Multiplies the current transformation on the right (!) by the given transformation.
    self_type operator*(const self_type& other) const;
    const self_type& operator*=(const self_type& other);

    // Resets current transformation to identity.
    void reset();

    // Returns the 4x4 matrix corresponding to the current transformation.
    TransformationMatrix matrix() const;

    // Allow stream operator<< access Transformation3D members.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const Transformation3D<V>& obj);

private:
    void set_rotation_block_(const Quaternion& q);
    void set_rotation_block_(const RotationMatrix& rot_matrix);
    void set_translation_block_(const Translation& t);

private:
    TransformationMatrix matrix_;
};

// Prints formatted transformation to the given stream.
template <typename RealType>
std::ostream& operator<<(std::ostream& os, const Transformation3D<RealType>& obj)
{
    // Print full vertex info.
    os << boost::format("Transformation3D object %1$#x, %2% bytes: ") % &obj % sizeof(obj)
       << std::endl << boost::format("    %1% %|22t|%2% %|40t|%3% %|58t|%4%")
          % obj.matrix_(0, 0) % obj.matrix_(0, 1) % obj.matrix_(0, 2) % obj.matrix_(0, 3)
       << std::endl << boost::format("    %1% %|22t|%2% %|40t|%3% %|58t|%4%")
          % obj.matrix_(1, 0) % obj.matrix_(1, 1) % obj.matrix_(1, 2) % obj.matrix_(1, 3)
       << std::endl << boost::format("    %1% %|22t|%2% %|40t|%3% %|58t|%4%")
          % obj.matrix_(2, 0) % obj.matrix_(2, 1) % obj.matrix_(2, 2) % obj.matrix_(2, 3)
       << std::endl << boost::format("    %1% %|22t|%2% %|40t|%3% %|58t|%4%")
          % obj.matrix_(3, 0) % obj.matrix_(3, 1) % obj.matrix_(3, 2) % obj.matrix_(3, 3)
       << std::endl << boost::format("end of object %1$#x.") % &obj << std::endl;

    return os;
}

template <typename RealType> inline
Transformation3D<RealType>::Transformation3D()
{
    reset();
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const Quaternion& q)
{
    reset();
    set_rotation_block_(q);
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const RotationMatrix& rot_matrix)
{
    reset();
    set_rotation_block_(rot_matrix);
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const Translation& t)
{
    reset();
    set_translation_block_(t); 
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const Quaternion& q, const Translation& t)
{
    reset();
    set_rotation_block_(q);
    set_translation_block_(t);
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const RotationMatrix& rot_matrix, const Translation& t)
{
    reset();
    set_rotation_block_(rot_matrix);
    set_translation_block_(t);
}

template <typename RealType>
Transformation3D<RealType>::Transformation3D(const TransformationMatrix& trans_matrix)
{
    reset();
    matrix_ = trans_matrix;
}


template <typename RealType> inline
typename Transformation3D<RealType>::Point3D Transformation3D<RealType>::operator*(
    Point3D point) const
{
    // Convert Point3D to boost BLAS vector. Convert to homogeneous coordinates.
    math::bounded_vector<RealType, 4> homog_vec;
    math::subrange(homog_vec, 0, 3) = math::from_bo_vector(point);

    // Perform multiplication.
    boost::numeric::ublas::bounded_vector<RealType, 4> result = math::prod(matrix_, homog_vec);

    // Convert back from boost BLAS vector to Point3D. Convert back from homogeneous
    // coordinates.
//    Point3D retvalue = math::to_bo_vector(math::subrange(result, 0, 3));
    math::bounded_vector<RealType, 3> sub(math::subrange(result, 0, 3));
    Point3D retvalue = math::to_bo_vector(sub);

    return
        retvalue;
}

template <typename RealType> inline
Transformation3D<RealType> Transformation3D<RealType>::operator*(const self_type& other) const
{
    self_type retvalue(*this);
    retvalue *= other;
    return retvalue;
}

template <typename RealType> inline
const Transformation3D<RealType>& Transformation3D<RealType>::operator*=(const self_type& other)
{
    matrix_ = math::prod(matrix_, other.matrix_);
    return *this;
}

template <typename RealType> inline
void Transformation3D<RealType>::reset()
{
    matrix_ = math::identity_matrix<RealType>(4);
}

template <typename RealType>
typename Transformation3D<RealType>::TransformationMatrix Transformation3D<RealType>::matrix() const
{
    return matrix_;
}

template <typename RealType>
void Transformation3D<RealType>::set_rotation_block_(const Quaternion& q)
{
    // Precompute some values based on quaternion components.
    RealType q00 = q[0]*q[0];
    RealType q11 = q[1]*q[1];
    RealType q22 = q[2]*q[2];
    RealType q33 = q[3]*q[3];
    RealType q03 = q[0]*q[3];
    RealType q13 = q[1]*q[3];
    RealType q23 = q[2]*q[3];
    RealType q02 = q[0]*q[2];
    RealType q12 = q[1]*q[2];
    RealType q01 = q[0]*q[1];

    // Fill the top-left 3x3 block (corresponding to rotation).
    matrix_(0, 0) = q00 + q11 - q22 - q33;
    matrix_(1, 1) = q00 - q11 + q22 - q33;
    matrix_(2, 2) = q00 - q11 - q22 + q33;
    matrix_(0, 1) = RealType(2) * (q12 - q03);
    matrix_(1, 0) = RealType(2) * (q12 + q03);
    matrix_(0, 2) = RealType(2) * (q13 + q02);
    matrix_(2, 0) = RealType(2) * (q13 - q02);
    matrix_(1, 2) = RealType(2) * (q23 - q01);
    matrix_(2, 1) = RealType(2) * (q23 + q01);
}

template <typename RealType>
void Transformation3D<RealType>::set_rotation_block_(const RotationMatrix& rot_matrix)
{
    // Fill the top left 3x3 block (corresponding to rotation).
    math::subrange(matrix_, 0, 3, 0, 3) = rot_matrix;
}

template <typename RealType>
void Transformation3D<RealType>::set_translation_block_(const Translation& t)
{
    // Fill the last column (corresponding to translation) using given translation vector.
//    math::subrange(math::column(matrix_, 3), 0, 3) = t;
    matrix_(0, 3) = t[0];
    matrix_(1, 3) = t[1];
    matrix_(2, 3) = t[2];
}

} // namespace math
} // namespace bo

#endif // TRANSFORMATION_3D_HPP_23877F3F_5EB9_4221_A422_C239EC9AA5D6_
