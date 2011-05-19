
#include "gtest/gtest.h"

#include "common/vector.hpp"
using namespace common;


// Create a so-called "text fixture" using base class form GTEST.
class VectorTest: public testing::Test
{
protected:

    VectorTest(): vec2_(5.f), vec3_(0.3, 4., 17), vec4_(9)
    { }

    // Variables for all Vector<> tests to use.
    Vector<float, 3> vec1_;
    Vector<double, 3> vec2_;
    Vector<double, 3> vec3_;
    Vector<int, 4> vec4_;
};