
/******************************************************************************

    vector.hpp, v 1.1.2 2011.04.12

    Multidimensional Vector (Point) class. 

    Copyright (c) 2010, 2011
    Alexander Rukletsov <alexander.rukletsov@ziti.uni-heidelberg.de>
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

#ifndef VECTOR_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_
#define VECTOR_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_

#include <cstddef>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <boost/format.hpp>
#include <boost/static_assert.hpp>


namespace common {

// Basic n-vector. Some methods, such as x(), y(), cross_product() are available
// only for certain specializations. A compiler error will be thrown if an 
// inappropriate call attempt is made. See notes on accessor functions below.
template <typename T, std::size_t N>
class Vector
    : boost::additive2< Vector<T, N>, T
    , boost::multiplicative2< Vector<T, N>, T
    , boost::equality_comparable1< Vector<T, N>
    , boost::additive1< Vector<T, N>
    > > > > 
{
public:
    // Each component is initialized either to 0 or to a given value.
    // If a c-array of type S is given, copy its elements.
	Vector();
    Vector(const T& value);
    template <typename S> explicit Vector(const S& data, std::size_t length);

    // These constructors add support for specializations with N = 2, 3, 4. See 
    // notes on accessor functions below.
    Vector(const T& x, const T& y);
    Vector(const T& x, const T& y, const T& z);
    Vector(const T& x, const T& y, const T& z, const T& w);

    // Some of standard operators. boost::operators library adds more.
    const Vector<T, N>& operator+=(const T& scalar);
    const Vector<T, N>& operator-=(const T& scalar);
    const Vector<T, N>& operator*=(const T& scalar);
    const Vector<T, N>& operator/=(const T& scalar);

    // In general won't work for floats. This is because not every real number can 
    // be represented by float/double/long double and therefore theoretically equal
    // numbers can differ. I.e. f^{-1}(f(x)) can differ from x. Fore more information
    // on this topic see
    //     http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    bool operator==(const Vector<T, N>& other) const;

    const Vector<T, N>& operator+=(const Vector<T, N>& other);
    const Vector<T, N>& operator-=(const Vector<T, N>& other);

    // Assignment and access operators. Range-check is done by boost::array via 
    // debug-only assertions. Use at() method for safer but less efficient version 
    // with exceptions.
    const T& operator[](std::size_t index) const;
    T& operator[](std::size_t index);

    // Assignment and access methods. Throws an exception in case of bad index.
    // Safer, but less efficient alternative of opeartor[].
    const T& at(std::size_t index) const;
    T& at(std::size_t index);


    // These special accessor functions are available only where appropriate.
    // That means, if you try to call x() on Vector<T, 10> or z() on Vector<T, 2>
    // you'll get an error. Unfortunately, not "... is not a member of ..." error,
    // but still clear enough to understand the reason. 
    //
    // The other solution is to have a basic, say, VectorImpl<T, N> class and 
    // derived classes Vector<T, N> and its specializations Vector<T, 2> and so on.
    // But this approach leads to a lot of copy-paste methods, e.g. all operators,
    // constructors and some other. The possible solution is to write a macros for 
    // automatic creation of these duplicated methods. An open question is performance
    // in the presence of inheritance and possible redundant object copying. That 
    // should be checked in case of this approach.
    //
    // Finally, it was decided to favor the first approach because of its simplicity
    // and not worse performance.

    const T& x() const;
    const T& y() const;
    const T& z() const;
    const T& w() const;

    T& x();
    T& y();
    T& z();
    T& w();
    
    // Cross product makes sense only in 3D and therefore is available only for 
    // Vector<T, 3>.
    Vector<T, N> cross_product(const Vector<T, N>& other) const;

    // See below for operator* (dot product), defined outside the class.
    // See below for operator<<, defined outside the class.

    // Simple aggregation functions. 
    T min() const;
    T max() const;
    T sum() const;
    T product() const;
    T avg() const;
	
    // These functions use std::min_element() and therefore return an index of the 
    // first found element (in case there are several equal min/max elements).
    std::size_t min_index() const;
    std::size_t max_index() const;

    // Compute euclidean norm of the vectors. If the return variable is given,
    // try to convert the result to retvar's type.
    double eucl_norm() const;
    template <typename RetType> void eucl_norm(RetType& retvar) const;

    // Note that for integral types normalize won't work. For this reason this
    // function is designed const and it returns a normalized double vector.
    // If the vector is a null-vector, normalization is sensless. However, some
    // software (e.g. Wolfram Mathematica 7.0) returns the null-vector as a result
    // of a normalization of the null-vector. In order not to inflate the code,
    // the function called on the null-vector will make 0/0 division and return
    // a vector made of NaNs. 
    Vector<double, N> normalized() const;

    // Size is always the same: N.
    std::size_t size() const;

    // Swap method (linear complexity).
    void swap(Vector<T, N>& other);

    // Fill each component with a given value.
    void fill(const T& value);

    // Set components' values from anything, that can be accessed by [].
    template <typename S> void assign(const S& data, std::size_t length);

protected:
    // Provide range check for a given index. Throws a std::out_of_range exception
    // in case of bad index.
    void check_range(std::size_t index) const;

protected:
    boost::array<T, N> components;
};


// Dot product operator for two Vector<T, N>.
template <typename T, std::size_t N> inline
T operator*(Vector<T, N> lhs, const Vector<T, N>& rhs) 
{ 
    T retvalue = 0;
    for (std::size_t i = 0; i < N; ++i)
        retvalue += lhs[i] * rhs[i];
    	
    return retvalue;
}

// Stream operator<< for printing Vector<T, N> contents.
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream &os, const Vector<T, N>& obj)
{
    os << boost::format("%1%-Vector, object %2$#x, %3% bytes: ") 
        % N % &obj % sizeof(obj) << std::endl << "    (";

    for (std::size_t i = 0; i < N-1; ++i)
        os << boost::format("%1%, %|4t|") % obj[i];

    // Print last element separately in order to avoid last comma and spaces.
    os << boost::format("%1%)") % obj[N-1] << std::endl 
       << boost::format("end of object %1$#x.") % &obj << std::endl;

    return os;
}


template <typename T, std::size_t N>
Vector<T, N>::Vector()
{ 
    fill(static_cast<T>(0));
}

template <typename T, std::size_t N>
Vector<T, N>::Vector(const T& value)
{ 
    fill(value);
}

template <typename T, std::size_t N> template <typename S>
Vector<T, N>::Vector(const S& data, std::size_t length)
{
    assign(data, length);
}

template <typename T, std::size_t N>
Vector<T, N>::Vector(const T& x, const T& y)
{
    BOOST_STATIC_ASSERT(N == 2);

    components[0] = x;
    components[1] = y;
}

template <typename T, std::size_t N>
Vector<T, N>::Vector(const T& x, const T& y, const T& z)
{
    BOOST_STATIC_ASSERT(N == 3);

    components[0] = x;
    components[1] = y;
    components[2] = z;
}

template <typename T, std::size_t N>
Vector<T, N>::Vector(const T& x, const T& y, const T& z, const T& w)
{
    BOOST_STATIC_ASSERT(N == 4);

    components[0] = x;
    components[1] = y;
    components[2] = z;
    components[3] = w;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator+=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] += scalar;

    return *this;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator-=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] -= scalar;

    return *this;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator*=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] *= scalar;

    return *this;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator/=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] /= scalar;

    return *this;
}

template <typename T, std::size_t N>
bool Vector<T, N>::operator==(const Vector<T, N>& other) const
{
    for (std::size_t i = 0; i < N; ++i)
        if (components[i] != other.components[i])
            return false;

    return true;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator+=(const Vector<T, N>& other)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] += other.components[i];

    return *this;
}

template <typename T, std::size_t N>
const Vector<T, N>& Vector<T, N>::operator-=(const Vector<T, N>& other)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] -= other.components[i];

    return *this;
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::operator[](std::size_t index) const
{
    return components[index];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::operator[](std::size_t index)
{
    return components[index];
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::at(std::size_t index) const
{
    check_range(index);
    return components[index];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::at(std::size_t index)
{
    check_range(index);
    return components[index];
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::x() const
{
    BOOST_STATIC_ASSERT(N >= 1 && N <= 4);
    return components[0];
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::y() const
{
    BOOST_STATIC_ASSERT(N >= 2 && N <= 4);
    return components[1];
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::z() const
{
    BOOST_STATIC_ASSERT(N >= 3 && N <= 4);
    return components[2];
}

template <typename T, std::size_t N> inline
const T& Vector<T, N>::w() const
{
    BOOST_STATIC_ASSERT(N == 4);
    return components[3];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::x()
{
    BOOST_STATIC_ASSERT(N >= 1 && N <= 4);
    return components[0];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::y()
{
    BOOST_STATIC_ASSERT(N >= 2 && N <= 4);
    return components[1];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::z()
{
    BOOST_STATIC_ASSERT(N >= 3 && N <= 4);
    return components[2];
}

template <typename T, std::size_t N> inline
T& Vector<T, N>::w()
{
    BOOST_STATIC_ASSERT(N == 4);
    return components[3];
}

template <typename T, std::size_t N>
Vector<T, N> Vector<T, N>::cross_product(const Vector<T, N>& other) const
{
    BOOST_STATIC_ASSERT(N == 3);

    return Vector<T, N>
        (y() * other.z() - z() * other.y(),
         z() * other.x() - x() * other.z(),
         x() * other.y() - y() * other.x());
}

template <typename T, std::size_t N>
T Vector<T, N>::min() const
{ 
    return 
        (*std::min_element(components.begin(), components.end()));
}

template <typename T, std::size_t N>
std::size_t Vector<T, N>::min_index() const
{
    return 
        std::distance(components.begin(), 
            std::min_element(components.begin(), components.end())); 
}

template <typename T, std::size_t N>
T Vector<T, N>::max() const
{ 
    return 
        (*std::max_element(components.begin(), components.end()));
}

template <typename T, std::size_t N>
std::size_t Vector<T, N>::max_index() const
{
    return 
        std::distance(components.begin(), 
            std::max_element(components.begin(), components.end()));
}

template <typename T, std::size_t N> 
T Vector<T, N>::sum() const 
{ 
    // A small optimization here: start with the second elem and pass first elem
    // as an initial value.
    return
        std::accumulate(components.begin() + 1, components.end(), components[0]);
}

template <typename T, std::size_t N>
T Vector<T, N>::product() const
{ 
    // A small optimization here: start with the second elem and pass first elem
    // as an initial value.
    return
        std::accumulate(components.begin() + 1, components.end(), components[0], 
            std::multiplies<T>());
}

template <typename T, std::size_t N>
T Vector<T, N>::avg() const
{ 
    return 
        (sum() / N); 
}

template <typename T, std::size_t N> 
double Vector<T, N>::eucl_norm() const
{
    return 
        sqrt(static_cast<double>((*this) * (*this)));
}

template <typename T, std::size_t N> template <typename RetType> 
void Vector<T, N>::eucl_norm(RetType& retvar) const
{
    retvar = static_cast<RetType>(sqrt(static_cast<double>((*this) * (*this))));
}

template <typename T, std::size_t N> 
Vector<double, N> Vector<T, N>::normalized() const
{
    double factor = 1.0 / eucl_norm();
    Vector<double, N> retvalue;

    for (std::size_t i = 0; i < N; ++i)
        retvalue[i] = static_cast<double>(components[i]) * factor;

    return retvalue;
}

template <typename T, std::size_t N> inline
std::size_t Vector<T, N>::size() const
{
    return N;
}

template <typename T, std::size_t N> 
void Vector<T, N>::swap(Vector<T, N>& other)
{
    components.swap(other.components);
}

template <typename T, std::size_t N> 
void Vector<T, N>::fill(const T& value)
{
    components.fill(value);
}

template <typename T, std::size_t N> template <typename S>
void Vector<T, N>::assign(const S& data, std::size_t length)
{
    // If a given array is smaller than N, copy everything and set 0 for other
    // components. If a given array is bigger than N, copy first N elements.
    std::size_t num = (length < N) ? length : N;
    for (std::size_t i = 0; i < num; ++i)
        components[i] = static_cast<T>(data[i]);

    for (std::size_t i = num; i < N; ++i)
        components[i] = static_cast<T>(0);
}

template <typename T, std::size_t N> 
void Vector<T, N>::check_range(std::size_t index) const
{
    if (index >= size()) 
    {
        std::out_of_range e((boost::format(
            "%1%-Vector's index \"%2%\" is out of range.") % N % index).str());

        throw e;
    }
}

} // namespace common

#endif // VECTOR_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_
