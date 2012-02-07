
/******************************************************************************

    extended_math.hpp, v 1.0.0 2011.09.28

    Extension of the standard <cmath> header.

    Copyright (c) 2011
    Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_
#define EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_

#include <cmath>
#include <boost/math/tr1.hpp>

namespace common {

template <typename T> inline
T square(const T& arg)
{
    return (arg * arg);
}

template <typename T> inline
T cube(const T& arg)
{
    return (arg * arg * arg);
}

} // namespace common

#endif // EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_