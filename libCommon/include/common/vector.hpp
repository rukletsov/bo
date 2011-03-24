
/******************************************************************************

    vector.hpp, v 1.0.2 2011.03.23

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
#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <boost/format.hpp>


namespace common {

// Basic n-vector class. Should not be used directly. Use Vector<T, N> instead.
template <typename T, std::size_t N>
class VectorImpl
    : boost::additive2< VectorImpl<T, N>, T
    , boost::multiplicative2< VectorImpl<T, N>, T
    , boost::equality_comparable1< VectorImpl<T, N>
    , boost::additive1< VectorImpl<T, N>
    > > > > 
{

public:
    // Each component is initialized either to 0 or to a given value.
    // If a c-array of type S is given, copy its elements.
	VectorImpl();
    VectorImpl(const T& value);
    template <typename S> explicit VectorImpl(const S& data, std::size_t length);

    // Some of standard operators. boost::operators adds more.
    VectorImpl<T, N> operator+=(const T& scalar);
    VectorImpl<T, N> operator-=(const T& scalar);
    VectorImpl<T, N> operator*=(const T& scalar);
    VectorImpl<T, N> operator/=(const T& scalar);

    bool operator==(const VectorImpl<T, N>& other) const;

    VectorImpl<T, N> operator+=(const VectorImpl<T, N>& other);
    VectorImpl<T, N> operator-=(const VectorImpl<T, N>& other);

    // Dot product operator cannot be created by boost since its return value is
    // T, not VectorImpl<T, N>. Therefore, create these operators manually. See below
    // for operator* implementation.
    T operator*=(const VectorImpl<T, N>& other);  

    // Assignment and access operators. Range-check is done by boost::array.
    const T& operator[](std::size_t index) const;
    T& operator[](std::size_t index);

    // See below for operator<<, defined outside the class.

    // Simple usual functions.
    T min() const;
    std::size_t min_index() const;
    T max() const;
    std::size_t max_index() const;
    T sum() const;
    T product() const;
    T avg() const;
	
    // Compute euclidean norm of the vectors. If the return variable is given,
    // try to convert the result to retvar's type.
    double eucl_norm() const;
    template <typename RetType> void eucl_norm(RetType& retvar) const;

    // Note that for integral types normalize won't work. For this reason this
    // function is designed const and it returns a normalized double vector.
    VectorImpl<double, N> normalized() const;

    // Size is always the same: N.
    std::size_t size() const;

    // Swap method (linear complexity).
    void swap(VectorImpl<T, N>& other);

    // Fill each component with a given value.
    void fill(const T& value);

    // Set components' values from anything, that can be accessed by [].
    template <typename S> void assign(const S& data, std::size_t length);

protected:
    boost::array<T, N> components;
};


template <typename T, std::size_t N>
class Vector: public VectorImpl<T, N>
    , private boost::additive2< Vector<T, N>, T
    , boost::multiplicative2< Vector<T, N>, T
    , boost::equality_comparable1< Vector<T, N>
    , boost::additive1< Vector<T, N>
    > > > > 
{

public:
    // See VectorImpl<> for details.
	Vector();
    Vector(const T& value);
    template <typename S> explicit Vector(const S& data, std::size_t length);
};


template <typename T>
class Vector<T, 2>: public VectorImpl<T, 2>
{

public:
    // See VectorImpl<> for details.
	Vector();
    Vector(const T& value);

    // This version of constructor exists only in this template specialization.
    Vector(const T& x, const T& y);

    // Additional assignment and access functions. Range-check is done by boost::array.
    const T& x() const;
    const T& y() const;
    T& x();
    T& y();
};


template <typename T>
class Vector<T, 3>: public VectorImpl<T, 3>
{

public:
    // See VectorImpl<> for details.
	Vector();
    Vector(const T& value);

    // This version of constructor exists only in this template specialization.
    Vector(const T& x, const T& y, const T& z);

    // Additional assignment and access functions. Range-check is done by boost::array.
    const T& x() const;
    const T& y() const;
    const T& z() const;
    T& x();
    T& y();
    T& z();
};



// Dot product operator for two VectorImpl<T, N>.
template <typename T, std::size_t N> inline
T operator*(VectorImpl<T, N> lhs, const VectorImpl<T, N>& rhs) 
{ 
    return lhs *= rhs; 
}

// Stream operator<< for printing VectorImpl<T, N> contents.
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream &os, const VectorImpl<T, N>& obj)
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
VectorImpl<T, N>::VectorImpl()
{ 
    components.fill(static_cast<T>(0));
}

template <typename T, std::size_t N>
VectorImpl<T, N>::VectorImpl(const T& value)
{ 
    fill(value);
}

template <typename T, std::size_t N> template <typename S>
VectorImpl<T, N>::VectorImpl(const S& data, std::size_t length)
{
    assign(data, length);
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator+=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] += scalar;

    return *this;
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator-=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] -= scalar;

    return *this;
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator*=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] *= scalar;

    return *this;
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator/=(const T& scalar)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] /= scalar;

    return *this;
}

template <typename T, std::size_t N>
bool VectorImpl<T, N>::operator==(const VectorImpl<T, N>& other) const
{
    for (std::size_t i = 0; i < N; ++i)
        if (components[i] != other.components[i])
            return false;

    return true;
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator+=(const VectorImpl<T, N>& other)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] += other.components[i];

    return *this;
}

template <typename T, std::size_t N>
VectorImpl<T, N> VectorImpl<T, N>::operator-=(const VectorImpl<T, N>& other)
{
    for (std::size_t i = 0; i < N; ++i)
        components[i] -= other.components[i];

    return *this;
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::operator*=(const VectorImpl<T, N>& other)
{
    T retvalue = 0;
    for (std::size_t i = 0; i < N; ++i)
        retvalue += components[i] * other.components[i];
    	
    return retvalue;
}

template <typename T, std::size_t N> inline
const T& VectorImpl<T, N>::operator[](std::size_t index) const
{
    return components[index];
}

template <typename T, std::size_t N> inline
T& VectorImpl<T, N>::operator[](std::size_t index)
{
    return components[index];
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::min() const
{ 
    T min_value = components[0];
    for (std::size_t i = 1; i < N; ++i)
	    if (components[i] < min_value)  
            min_value = components[i];

    return min_value; 
}

template <typename T, std::size_t N>
std::size_t VectorImpl<T, N>::min_index() const
{
    std::size_t min_index = 0;
    T min_value = components[0];

    for (std::size_t i = 1; i < N; ++i)
    {
	    if (components[i] < min_value)  
        {
            min_index = i;
            min_value = components[i];
        }
    }

    return min_index;
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::max() const
{ 
    T max_value = components[0];
    for (std::size_t i = 1; i < N; ++i)
	    if (components[i] > max_value)  
            max_value = components[i];

    return max_value; 
}

template <typename T, std::size_t N>
std::size_t VectorImpl<T, N>::max_index() const
{
    std::size_t max_index = 0;
    T max_value = components[0];

    for (std::size_t i = 1; i < N; ++i)
    {
	    if (components[i] > max_value)  
        {
            max_index = i;
            max_value = components[i];
        }
    }

    return max_index;
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::sum() const
{ 
    T total = components[0];
    for (std::size_t i = 1; i < N; ++i) 
        total += components[i];

    return total; 
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::product() const
{ 
    T product = components[0];
    for (std::size_t i = 1; i < N; ++i) 
        product *= components[i];

    return product; 
}

template <typename T, std::size_t N>
T VectorImpl<T, N>::avg() const
{ 
    return 
        (sum() / N); 
}

template <typename T, std::size_t N> 
double VectorImpl<T, N>::eucl_norm() const
{
    return 
        sqrt(static_cast<double>((*this) * (*this)));
}

template <typename T, std::size_t N> template <typename RetType> 
void VectorImpl<T, N>::eucl_norm(RetType& retvar) const
{
    retvar = static_cast<RetType>(sqrt(static_cast<double>((*this) * (*this))));
}

template <typename T, std::size_t N> 
VectorImpl<double, N> VectorImpl<T, N>::normalized() const
{
    double factor = 1.0 / eucl_norm();
    VectorImpl<double, N> retvalue;

    for (std::size_t i = 0; i < N; ++i)
        retvalue[i] = static_cast<double>(components[i]) * factor;

    return retvalue;
}

template <typename T, std::size_t N> inline
std::size_t VectorImpl<T, N>::size() const
{
    return N;
}

template <typename T, std::size_t N> 
void VectorImpl<T, N>::swap(VectorImpl<T, N>& other)
{
    components.swap(other.components);
}

template <typename T, std::size_t N> 
void VectorImpl<T, N>::fill(const T& value)
{
    components.fill(value);
}

template <typename T, std::size_t N> template <typename S>
void VectorImpl<T, N>::assign(const S& data, std::size_t length)
{
    // If a given array is smaller than N, copy everything and set 0 for other
    // components. If a given array is bigger than N, copy first N elements.
    std::size_t num = length < N ? length : N;
    for (std::size_t i = 0; i < num; ++i)
        components[i] = static_cast<T>(data[i]);

    for (std::size_t i = num; i < N; ++i)
        components[i] = static_cast<T>(0);
}


// Vector<T, N> methods.
template <typename T, std::size_t N>
Vector<T, N>::Vector(): VectorImpl<T, N>()
{ }

template <typename T, std::size_t N>
Vector<T, N>::Vector(const T& value): VectorImpl<T, N>(value)
{ }

template <typename T, std::size_t N> template <typename S>
Vector<T, N>::Vector(const S& data, std::size_t length): VectorImpl<T, N>(data, length)
{ }


// Vector<T, 2> methods.
template <typename T>
Vector<T, 2>::Vector(): VectorImpl<T, 2>()
{ }

template <typename T>
Vector<T, 2>::Vector(const T& value): VectorImpl<T, 2>(value)
{ }

template <typename T> 
Vector<T, 2>::Vector(const T& x, const T& y)
{
    this->components[0] = x;
    this->components[1] = y;
}

template <typename T> inline
const T& Vector<T, 2>::x() const
{
    return components[0];
}

template <typename T> inline
const T& Vector<T, 2>::y() const
{
    return components[1];
}

template <typename T> inline
T& Vector<T, 2>::x() 
{
    return components[0];
}

template <typename T> inline
T& Vector<T, 2>::y() 
{
    return components[1];
}


// Vector<T, 3> methods.
template <typename T>
Vector<T, 3>::Vector(): VectorImpl<T, 3>()
{ }

template <typename T>
Vector<T, 3>::Vector(const T& value): VectorImpl<T, 3>(value)
{ }

template <typename T> 
Vector<T, 3>::Vector(const T& x, const T& y, const T& z)
{
    this->components[0] = x;
    this->components[1] = y;
    this->components[2] = z;
}

template <typename T> inline
const T& Vector<T, 3>::x() const
{
    return components[0];
}

template <typename T> inline
const T& Vector<T, 3>::y() const
{
    return components[1];
}

template <typename T> inline
const T& Vector<T, 3>::z() const
{
    return components[2];
}

template <typename T> inline
T& Vector<T, 3>::x() 
{
    return components[0];
}

template <typename T> inline
T& Vector<T, 3>::y() 
{
    return components[1];
}

template <typename T> inline
T& Vector<T, 3>::z() 
{
    return components[2];
}



























// 2D vector.
template<typename T>
struct Vector2
{
	T x;
	T y;

	Vector2();
	Vector2(T _x, T _y);
	
	T get_dimension_value(int dimension_index);
	void set_dimension_value(int dimension_index, T value);
};


// 3D vector.
template<typename T>
struct Vector3
{
	T x;
	T y;
	T z;

	Vector3();
	Vector3(T _x, T _y, T _z);

	Vector3<T> operator + (const Vector3<T>& other) const;
	Vector3<T> operator - (const Vector3<T>& other) const;
	Vector3<T> operator * (const T& scalar) const;
	Vector3<T> operator / (const T& scalar) const;
	T operator * (const Vector3<T>& other) const;
	bool operator == (const Vector3<T>& other) const;

    T dot_product(const Vector3<T>& other) const;
    Vector3<T> cross_product(const Vector3& other) const;

	double get_eucl_norm() const;

    // Note that for integral types normalize won't work. For this reason this
    // function is designed const and it returns a normalized double vector.
    Vector3<double> normalized() const;

	T get_dimension_value(int dimension_index) const;
	void set_dimension_value(int dimension_index, T value);
};

 
// 4D vector.
template<typename T>
struct Vector4
{
	T x1;
	T y1;
	T x2;
	T y2;

	Vector4();
	Vector4(T _x1, T _y1, T _x2, T _y2);
	Vector4(const Vector2<T> &p1, const Vector2<T> &p2);

	T get_dimension_value(int dimension_index);
	void set_dimension_value(int dimension_index, T value);
};


template<typename T>
Vector2<T>::Vector2(): x(0), y(0)
{ } 

template<typename T>
Vector2<T>::Vector2(T _x, T _y): x(_x), y(_y)
{ }

template<typename T>
T Vector2<T>::get_dimension_value(int dimension_index)
{
	T dimension_value;
	switch (dimension_index)
	{
		case 0:
			dimension_value = this->x;
			break;
		case 1:
			dimension_value = this->y;
			break;
	}

	return dimension_value;
}

template<typename T>
void Vector2<T>::set_dimension_value(int dimension_index, T value)
{
	switch (dimension_index)
	{
		case 0:
			this->x=value;	
			break;
		case 1:
			this->y=value;
			break;
	}
}


template<typename T>
Vector3<T>::Vector3(): x(0), y(0), z(0)
{ }

template<typename T>
Vector3<T>::Vector3(T _x, T _y, T _z): x(_x), y(_y), z(_z)
{ }

template<typename T>
T Vector3<T>::get_dimension_value(int dimension_index) const
{
	T dimension_value;
	switch (dimension_index)
	{
		case 0:
			dimension_value = this->x;
			break;
		case 1:
			dimension_value = this->y;
			break;
		case 2:
			dimension_value = this->z;
            break;
	}

	return dimension_value;
}

template<typename T>
void Vector3<T>::set_dimension_value(int dimension_index, T value)
{
	switch (dimension_index)
	{
		case 0:
			this->x=value;	
			break;
		case 1:
			this->y=value;
			break;
		case 2:
			this->z=value;
			break;
	}
}

template<typename T>
Vector3<T> Vector3<T>::operator+(const Vector3<T>& other) const
{
	Vector3<T> p;
	p.x = this->x + other.x;
	p.y = this->y + other.y;
	p.z = this->z + other.z;

	return p;
}

template<typename T>
Vector3<T> Vector3<T>::operator-(const Vector3<T>& other) const
{
	Vector3<T> p;
	p.x = this->x - other.x;
	p.y = this->y - other.y;
	p.z = this->z - other.z;

	return p;
}

template<typename T>
Vector3<T> Vector3<T>::operator*(const T& scalar) const
{
	Vector3<T> p;
	p.x = this->x * scalar;
	p.y = this->y * scalar;
	p.z = this->z * scalar;

	return p;
}

template<typename T>
Vector3<T> Vector3<T>::operator/(const T& scalar) const
{
	Vector3<T> p;
	p.x = this->x / scalar;
	p.y = this->y / scalar;
	p.z = this->z / scalar;

	return p;
}

template<typename T>
T Vector3<T>::operator*(const Vector3<T>& other) const
{
	return 
        x * other.x + y * other.y + z * other.z;
}

template<typename T>
bool Vector3<T>::operator==(const Vector3<T>& other) const
{
    return 
        (((x == other.x) && (y==other.y) && (z == other.z)) ? true : false);
}

template<typename T>
T Vector3<T>::dot_product(const Vector3<T>& other) const
{
	return
        x * other.x + y * other.y + z * other.z;
}

template<typename T>
Vector3<T> Vector3<T>::cross_product(const Vector3& other) const
{
    return Vector3<T>
        (y * other.z - z * other.y,
         z * other.x - x * other.z,
         x * other.y - y * other.x);
}

template<typename T> 
double Vector3<T>::get_eucl_norm() const
{
    return 
        sqrt(static_cast<double>((*this) * (*this)));
}


template<typename T> 
Vector3<double> Vector3<T>::normalized() const
{
    double factor = 1.0 / get_eucl_norm();
    return Vector3<double>
        (static_cast<double>(x) * factor,
         static_cast<double>(y) * factor,
         static_cast<double>(z) * factor);
}


template<typename T>
Vector4<T>::Vector4(): x1(0), y1(0), x2(0), y2(0)
{ }

template<typename T>
Vector4<T>::Vector4(T _x1, T _y1, T _x2, T _y2): x1(_x1), y1(_y1), x2(_x2), y2(_y2)
{ }

template<typename T>
Vector4<T>::Vector4(const Vector2<T>& p1, const Vector2<T>& p2): x1(p1.x), y1(p1.y), 
    x2(p2.x), y2(p2.y)
{ }

template<typename T>
T Vector4<T>::get_dimension_value(int dimension_index)
{
	T dimension_value;
	switch (dimension_index)
	{
		case 0:
			dimension_value = this->x1;
			break;
		case 1:
			dimension_value = this->y1;
			break;
		case 2:
			dimension_value = this->x2;
			break;
		case 3:
			dimension_value = this->y2;
			break;
	}

	return dimension_value;
}

template<typename T>
void Vector4<T>::set_dimension_value(int dimension_index, T value)
{
	switch (dimension_index)
	{
		case 0:
			this->x1=value;	
			break;
		case 1:
			this->y1=value;
			break;
		case 2:
			this->x2=value;	
			break;
		case 3:
			this->y2=value;
	}
}

} // namespace common

#endif // VECTOR_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_
