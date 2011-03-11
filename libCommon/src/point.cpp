
#include "stdafx.h"

#include <math.h>

#include "point.hpp"


namespace common {

template<typename T>
Point2<T>::Point2(): x(0), y(0)
{ } 

template<typename T>
Point2<T>::Point2(T _x, T _y): x(_x), y(_y)
{ }

template<typename T>
T Point2<T>::get_dimension_value(int dimension_index)
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
void Point2<T>::set_dimension_value(int dimension_index, T value)
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
Point3<T>::Point3(): x(0), y(0), z(0)
{ }

template<typename T>
Point3<T>::Point3(T _x, T _y, T _z): x(_x), y(_y), z(_z)
{ }

template<typename T>
T Point3<T>::get_dimension_value(int dimension_index)
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
void Point3<T>::set_dimension_value(int dimension_index, T value)
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
Point3<T> Point3<T>::operator+(const Point3<T> &other) const
{
	Point3<T> p;
	p.x=this->x+other.x;
	p.y=this->y+other.y;
	p.z=this->z+other.z;

	return p;
}

template<typename T>
Point3<T> Point3<T>::operator-(const Point3<T> &other) const
{
	Point3<T> p;
	p.x=this->x-other.x;
	p.y=this->y-other.y;
	p.z=this->z-other.z;

	return p;
}

template<typename T>
Point3<T> Point3<T>::operator*(const T &scalar) const
{
	Point3<T> p;
	p.x=this->x*scalar;
	p.y=this->y*scalar;
	p.z=this->z*scalar;

	return p;
}

template<typename T>
Point3<T> Point3<T>::operator/(const T &scalar) const
{
	Point3<T> p;
	p.x=this->x/scalar;
	p.y=this->y/scalar;
	p.z=this->z/scalar;

	return p;
}

template<typename T>
T Point3<T>::operator*(const Point3<T> &other) const
{
	return x * other.x + y * other.y + z * other.z;
}

template<typename T>
bool Point3<T>::operator==(const Point3<T> &other) const
{
    return 
        (((x == other.x) && (y==other.y) && (z == other.z)) ? true : false);
}

template<typename T> 
double Point3<T>::get_eucl_norm() const
{
    return sqrt(double(x * x + y * y + z * z));
}


template<typename T>
Point4<T>::Point4(): x1(0), y1(0), x2(0), y2(0)
{ }

template<typename T>
Point4<T>::Point4(T _x1, T _y1, T _x2, T _y2): x1(_x1), y1(_y1), x2(_x2), y2(_y2)
{ }

template<typename T>
Point4<T>::Point4(const Point2<T> &p1, const Point2<T> &p2): x1(p1.x), y1(p1.y), 
    x2(p2.x), y2(p2.y)
{ }

template<typename T>
T Point4<T>::get_dimension_value(int dimension_index)
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
void Point4<T>::set_dimension_value(int dimension_index, T value)
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

template struct Point2<char>;
template struct Point2<unsigned char>;
template struct Point2<short>;
template struct Point2<unsigned short>;
template struct Point2<int>;
template struct Point2<unsigned int>;
template struct Point2<float>;
template struct Point2<double>;

template struct Point3<char>;
template struct Point3<unsigned char>;
template struct Point3<short>;
template struct Point3<unsigned short>;
template struct Point3<int>;
template struct Point3<unsigned int>;
template struct Point3<float>;
template struct Point3<double>;

template struct Point4<char>;
template struct Point4<unsigned char>;
template struct Point4<short>;
template struct Point4<unsigned short>;
template struct Point4<int>;
template struct Point4<unsigned int>;
template struct Point4<float>;
template struct Point4<double>;

} // namespace common
