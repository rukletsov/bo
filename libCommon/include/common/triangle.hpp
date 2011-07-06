
/******************************************************************************

    triangle.hpp, v 1.0.3 2011.03.22

    Triangle class. 

    Copyright (c) 2010, 2011
    Alexander Rukletsov <rukletsov@gmail.com>
    Dzmitry Hlindzich <dzmitry.hlindzich@ziti.uni-heidelberg.de>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
	    notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
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

#ifndef TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_
#define TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_

#include <stdexcept>
#include <boost/array.hpp>
#include <boost/operators.hpp>


namespace common {

template <typename PointType>
class Triangle: boost::equality_comparable1< Triangle<PointType> >
{ 

public:
    Triangle() { }
    Triangle(const PointType& _A, const PointType& _B, const PointType& _C);

    PointType A() const;
    PointType B() const;
    PointType C() const;

    // Assignment and access operators. Range-check is done by boost::array via
    // debug-only assertions. Use at() method for safer but less efficient version
    // with exceptions.
    const PointType& operator[](std::size_t index) const;
    PointType& operator[](std::size_t index);

    // Assignment and access methods. Throw an exception in case of bad index.
    // Safer, but less efficient alternative of opeartor[].
    const PointType& at(std::size_t index) const;
    PointType& at(std::size_t index);

    // In general won't work for floats. This is because not every real number can
    // be represented by float/double/long double and therefore theoretically equal
    // numbers can differ, i.e. f^{-1}(f(x)) can differ from x. Fore more information
    // on this topic see
    //     http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    bool operator==(const Triangle<PointType>& other) const;

protected:
    // Provides a range check for a given index. Throws a std::out_of_range exception
    // in case of bad index.
    static void check_range(std::size_t index) ;

protected:
    boost::array<PointType, 3> vertices_;
};


template <typename PointType>
Triangle<PointType>::Triangle(const PointType& _A, const PointType& _B,
                              const PointType& _C)
{
    vertices_[0] = _A;
    vertices_[1] = _B;
    vertices_[2] = _C;
}

template <typename PointType> inline
PointType Triangle<PointType>::A() const
{
    return vertices_[0];
}

template <typename PointType> inline
PointType Triangle<PointType>::B() const
{
    return vertices_[1];
}

template <typename PointType> inline
PointType Triangle<PointType>::C() const
{
    return vertices_[2];
}

template <typename PointType> inline
const PointType& Triangle<PointType>::operator[](std::size_t index) const
{
    return vertices_[index];
}

template <typename PointType> inline
PointType& Triangle<PointType>::operator[](std::size_t index)
{
    return vertices_[index];
}

template <typename PointType> inline
const PointType& Triangle<PointType>::at(std::size_t index) const
{
    check_range(index);
    return vertices_[index];
}

template <typename PointType> inline
PointType& Triangle<PointType>::at(std::size_t index)
{
    check_range(index);
    return vertices_[index];
}

template <typename PointType>
bool Triangle<PointType>::operator==(const Triangle<PointType>& other) const
{
    bool equal = !((vertices_[0] != other.vertices_[0]) ||
                   (vertices_[1] != other.vertices_[1]) ||
                   (vertices_[2] != other.vertices_[2]));
    return equal;
}

template <typename PointType>
void Triangle<PointType>::check_range(std::size_t index)
{
    if (index >= 3)
        throw std::out_of_range("Triangle has only 3 vertices.");
}

} // namespace common

#endif // TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_
