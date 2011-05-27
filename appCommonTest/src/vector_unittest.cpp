
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
    Vector<int, 3> int_vec(arr1, 2);

    EXPECT_EQ(int(1), int_vec[0]);
    EXPECT_EQ(int(2), int_vec[1]);
    EXPECT_EQ(int(0), int_vec[2]);

    int arr2[5] = {3, 4, 10, 45};
    Vector<double, 3> double_vec(arr2, 5);

    EXPECT_EQ(3., double_vec[0]);
    EXPECT_EQ(4., double_vec[1]);
    EXPECT_EQ(10., double_vec[2]);
}

TEST_F(VectorTest, CopyConstructor)
{
    Vector<double, 3> double_vec = vec2_;
    int diff = memcmp(&double_vec, &vec2_, 3 * sizeof(double));

    EXPECT_EQ(0, diff);
}

TEST_F(VectorTest, MemoryConsumption)
{
    // Vector is expected to consume no additional memory.
    EXPECT_EQ(3 * sizeof(float), sizeof(vec1_));
    EXPECT_EQ(4 * sizeof(int), sizeof(vec4_));
    EXPECT_EQ(sizeof(char), sizeof(Vector<int, 0>));
}

TEST_F(VectorTest, BoundaryChecks)
{
    EXPECT_NO_THROW(vec1_.at(0));
    EXPECT_NO_THROW(vec1_.at(vec1_.size() - 1));
    EXPECT_THROW(vec1_.at(vec1_.size()), std::out_of_range);
    EXPECT_THROW(vec1_.at(-1), std::out_of_range);

    // operator[] cannot be checked since it uses BOOST_ASSERT.
}

TEST_F(VectorTest, VectorScalarArithmetics)
{
    Vector<int, 4> int_vec = vec4_ * int(2);
    EXPECT_EQ(18, int_vec.x());
    EXPECT_EQ(18, int_vec.y());
    EXPECT_EQ(18, int_vec.z());
    EXPECT_EQ(18, int_vec.w());

    int_vec += int(2);
    EXPECT_EQ(20, int_vec.x());
    EXPECT_EQ(20, int_vec.y());
    EXPECT_EQ(20, int_vec.z());
    EXPECT_EQ(20, int_vec.w());

    int_vec = vec4_ - int(4);
    EXPECT_EQ(5, int_vec.x());
    EXPECT_EQ(5, int_vec.y());
    EXPECT_EQ(5, int_vec.z());
    EXPECT_EQ(5, int_vec.w());

    int_vec = vec4_ / int(4);
    EXPECT_EQ(2, int_vec.x());
    EXPECT_EQ(2, int_vec.y());
    EXPECT_EQ(2, int_vec.z());
    EXPECT_EQ(2, int_vec.w());
}

TEST_F(VectorTest, VectorVectorArithmetics)
{
    Vector<double, 3> double_vec = vec2_ + vec3_;
    EXPECT_DOUBLE_EQ(5.3, double_vec.x());
    EXPECT_DOUBLE_EQ(9., double_vec.y());
    EXPECT_DOUBLE_EQ(22., double_vec.z());

    double_vec -= vec2_;
    EXPECT_DOUBLE_EQ(vec3_.x(), double_vec.x());
    EXPECT_DOUBLE_EQ(vec3_.y(), double_vec.y());
    EXPECT_DOUBLE_EQ(vec3_.z(), double_vec.z());

    int arr1[] = {7, 8, 9};
    Vector<int, 4> int_vec(arr1, 3);
    int_vec.x() = int_vec.y() * 2 - int_vec.x();
    int_vec.y() = int_vec.y() * 9 / 8;
    int_vec.w() = int_vec.z();
    EXPECT_EQ(vec4_, int_vec);
}

TEST_F(VectorTest, DotProduct)
{
    double dot_product = vec2_ * vec3_;
    EXPECT_DOUBLE_EQ(106.5, dot_product);
}

TEST_F(VectorTest, CrossProduct)
{
    vec2_ = vec2_.cross_product(vec3_);
    EXPECT_DOUBLE_EQ(65., vec2_.x());
    EXPECT_DOUBLE_EQ(-83.5, vec2_.y());
    EXPECT_DOUBLE_EQ(18.5, vec2_.z());
}
