
/******************************************************************************

    triangle.hpp, v 1.0.2 2011.03.10

    Triangle class. 

    Copyright (c) 2010, 2011
    Alexander Rukletsov <alexander.rukletsov@ziti.uni-heidelberg.de>
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

#include <boost/array.hpp>


namespace common {

template <typename PointType>
class Triangle
{ 

public:
	Triangle() { };
	Triangle(const PointType& _A, const PointType& _B, const PointType& _C);

    PointType A() const;
	PointType B() const;
	PointType C() const;

    PointType& operator[] (std::size_t index);
    const PointType& operator[] (std::size_t index) const;

private:
    boost::array<PointType, 3> vertices;
};


template <typename PointType>
Triangle<PointType>::Triangle(const PointType& _A, const PointType& _B, const PointType& _C)
{
    vertices[0] = _A;
    vertices[1] = _B;
    vertices[2] = _C;
}

template <typename PointType> inline
PointType Triangle<PointType>::A() const
{
    return vertices[0];
}

template <typename PointType> inline
PointType Triangle<PointType>::B() const
{
    return vertices[1];
}

template <typename PointType> inline
PointType Triangle<PointType>::C() const
{
    return vertices[2];
}

template <typename PointType> inline
PointType& Triangle<PointType>::operator [](std::size_t index)
{
    return vertices[index];
}

template <typename PointType> inline
const PointType& Triangle<PointType>::operator [](std::size_t index) const
{
    return vertices[index];
}

} // namespace common

#endif // TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_
