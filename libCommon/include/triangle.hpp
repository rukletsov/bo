
#ifndef TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_
#define TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_

#include <boost/array.hpp>

#include "point.hpp"


namespace common {

template <typename T>
class Triangle2
{ 

public:
	Triangle2()
    { };
	Triangle2(const Point2<T>& _A, const Point2<T>& _B, const Point2<T>& _C);

    Point2<T> A() const;
	Point2<T> B() const;
	Point2<T> C() const;

    T& operator[] (std::size_t index);
    const T& operator[] (std::size_t index) const;

private:
    boost::array<Point2<T>, 3> vertices;
};

template <typename T>
class Triangle3
{

public:
	Triangle3()
    { };
	Triangle3(const Point3<T>& _A, const Point3<T>& _B, const Point3<T>& _C);

    Point3<T> A() const;
	Point3<T> B() const;
	Point3<T> C() const;

    T& operator[] (std::size_t index);
    const T& operator[] (std::size_t index) const;

private:
    boost::array<Point3<T>, 3> vertices;
};


template <typename T>
Triangle2<T>::Triangle2(const Point2<T>& _A, const Point2<T>& _B, const Point2<T>& _C)
{
    vertices[0] = _A;
    vertices[1] = _B;
    vertices[0] = _C;
}

template <typename T> inline
Point2<T> Triangle2<T>::A() const
{
    return vertices[0];
}

template <typename T> inline
Point2<T> Triangle2<T>::B() const
{
    return vertices[1];
}

template <typename T> inline
Point2<T> Triangle2<T>::C() const
{
    return vertices[2];
}

template <typename T> inline
T& Triangle2<T>::operator [](std::size_t index)
{
    return vertices[index];
}

template <typename T> inline
const T& Triangle2<T>::operator [](std::size_t index) const
{
    return vertices[index];
}


template <typename T>
Triangle3<T>::Triangle3(const Point3<T>& _A, const Point3<T>& _B, const Point3<T>& _C)
{
    vertices[0] = _A;
    vertices[1] = _B;
    vertices[0] = _C;
}

template <typename T> inline
Point3<T> Triangle3<T>::A() const
{
    return vertices[0];
}

template <typename T> inline
Point3<T> Triangle3<T>::B() const
{
    return vertices[1];
}

template <typename T> inline
Point3<T> Triangle3<T>::C() const
{
    return vertices[2];
}

template <typename T> inline
T& Triangle3<T>::operator [](std::size_t index)
{
    return vertices[index];
}

template <typename T> inline
const T& Triangle3<T>::operator [](std::size_t index) const
{
    return vertices[index];
}

} // namespace common

#endif // TRIANGLE_HPP_507AFC96_F3F4_40FF_827C_66F388AEDAD2_
