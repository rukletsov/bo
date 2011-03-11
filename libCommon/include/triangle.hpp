
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
    vertices[0] = _C;
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
