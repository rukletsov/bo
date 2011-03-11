
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
