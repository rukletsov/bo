
/******************************************************************************

    point.hpp, v 1.0.1 2011.03.10

    Multidimensional Point classes. 

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

#ifndef POINT_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_
#define POINT_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_


namespace common {

// 2D point.
template<typename T>
struct Point2
{
	T x;
	T y;

	Point2();
	Point2(T _x, T _y);
	
	T get_dimension_value(int dimension_index);
	void set_dimension_value(int dimension_index, T value);
};


// 3D point.
template<typename T>
struct Point3
{
	T x;
	T y;
	T z;

	Point3();
	Point3(T _x, T _y, T _z);

	Point3<T> operator + (const Point3<T> &other) const;
	Point3<T> operator - (const Point3<T> &other) const;
	Point3<T> operator * (const T &scalar) const;
	Point3<T> operator / (const T &scalar) const;
	T operator * (const Point3<T> &other) const;
	bool operator == (const Point3<T> &other) const;

	double get_eucl_norm() const;

	T get_dimension_value(int dimension_index);
	void set_dimension_value(int dimension_index, T value);
};


// 4D point.
template<typename T>
struct Point4
{
	T x1;
	T y1;
	T x2;
	T y2;

	Point4();
	Point4(T _x1, T _y1, T _x2, T _y2);
	Point4(const Point2<T> &p1, const Point2<T> &p2);

	T get_dimension_value(int dimension_index);
	void set_dimension_value(int dimension_index, T value);
};

} // namespace common

#endif // POINT_HPP_4545D406_43E3_4444_8E4B_9B5A10E7AB16_
