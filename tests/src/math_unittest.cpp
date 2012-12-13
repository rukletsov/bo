
#include <vector>
#include <gtest/gtest.h>

#include "bo/vector.hpp"
#include "bo/extended_math.hpp"
#include "debug_alloc.hpp"

using namespace bo;


// Create a so-called "text fixture" using base class form GTEST.
class MathTest: public testing::Test
{
public:
    typedef Vector<float, 3> Vec3f;
    typedef Vector<double, 5> Vec5d;
    typedef Vector<int, 2> Vec2i;

    typedef std::vector<Vec3f> Vecs3f;
    typedef std::vector<Vec5d> Vecs5d;
    typedef std::vector<Vec2i> Vecs2i;
    typedef std::vector<float> Samples1f;
    typedef std::vector<double> Samples1d;

protected:

    MathTest(): data3f_(4), data1d_(4)
    { }

    virtual void SetUp()
    {
        data3f_[0] = Vec3f(1.f, 2.f, 3.f); data3f_[1] = Vec3f(-1.f, -2.f, -3.f);
        data3f_[2] = Vec3f(2.f, 1.f, -3.f); data3f_[3] = Vec3f(-2.f, -1.f, 3.f);

        data1d_[0] = 1.2; data1d_[1] = 2.3; data1d_[2] = 3.4; data1d_[3] = 4.5;
    }

    // Data samples for tests.
    Vecs3f data3f_;
    Vecs5d data5d_;
    Vecs2i data2i_;
    Samples1f data1f_;
    Samples1d data1d_;
};


TEST_F(MathTest, Mean)
{
    Vec3f res3f = mean(data3f_);
    EXPECT_DOUBLE_EQ(0.f, res3f[0]);
    EXPECT_DOUBLE_EQ(0.f, res3f[1]);
    EXPECT_DOUBLE_EQ(0.f, res3f[2]);

    double res1d = mean(data1d_);
    EXPECT_DOUBLE_EQ(2.85, res1d);
}
