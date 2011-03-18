
/******************************************************************************

    vector.hpp, v 1.0.1 2011.03.10

    Multidimensional Vector and Point classes. 

    Copyright (c) 2010, 2011
    Dzmitry Hlindzich <dzmitry.hlindzich@ziti.uni-heidelberg.de>
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


namespace common {

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

	T get_dimension_value(int dimension_index);
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
T Vector3<T>::get_dimension_value(int dimension_index)
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
    return sqrt(double(x * x + y * y + z * z));
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
