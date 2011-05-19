
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

TEST_F(VectorTest, DefaultConstructor)
{
    // Default constructor should fill up the vector with zeros.

    float sum = 0.f;
    for (size_t i = 0; i < 3; ++i)
        sum += vec1_[i];

    EXPECT_EQ(0.f, sum);

    Vector<int, 3> zero_vec;
    int zeros[3];
    memset(&zeros, 0, 3 * sizeof(int));
    int diff = memcmp(&zero_vec[0], &zeros, 3 * sizeof(int));

    EXPECT_EQ(0, diff);
}

TEST_F(VectorTest, FillConstructor)
{
    double sum = 0.;
    for (size_t i = 0; i < 3; ++i)
        sum += vec2_[i];

    EXPECT_EQ(15., sum);
}

TEST_F(VectorTest, ArrayConstructor)
{
    float arr1[2] = {1.f, 2.2f};
    common::Vector<int, 3> int_vec(arr1, 2);

    EXPECT_EQ(int(1), int_vec[0]);
    EXPECT_EQ(int(2), int_vec[1]);
    EXPECT_EQ(int(0), int_vec[2]);

    int arr2[5] = {3, 4, 10, 45};
    common::Vector<double, 3> double_vec(arr2, 5);

    EXPECT_EQ(3., double_vec[0]);
    EXPECT_EQ(4., double_vec[1]);
    EXPECT_EQ(10., double_vec[2]);
}

TEST_F(VectorTest, CopyConstructor)
{
}

TEST_F(VectorTest, MemoryConsumption)
{
    // Vector should consume no additional memory.
    EXPECT_EQ(3 * sizeof(float), sizeof(vec1_));
    EXPECT_EQ(4 * sizeof(int), sizeof(vec4_));
    EXPECT_EQ(sizeof(char), sizeof(Vector<int, 0>));
}

