
#include <vector>
#include <gtest/gtest.h>

#include "bo/vector.hpp"
#include "bo/extended_math.hpp"
#include "bo/blas/pca.hpp"
#include "debug_alloc.hpp"

using namespace bo;

// Create a so-called "text fixture" using base class form GTEST.
class MathTest: public testing::Test
{
public:
    typedef Vector<float, 3> Vec3f;
    typedef Vector<double, 2> Vec2d;
    typedef Vector<int, 2> Vec2i;

    typedef std::vector<Vec3f> Vecs3f;
    typedef std::vector<Vec2d> Vecs2d;
    typedef std::vector<Vec2i> Vecs2i;
    typedef std::vector<float> Samples1f;
    typedef std::vector<double> Samples1d;

    typedef std::logic_error ex_t;

protected:

    MathTest(): data3f_(4), data2d_(10), data1d_(4)
    { }

    virtual void SetUp()
    {
        data3f_[0] = Vec3f(1.f, 2.f, 3.f); data3f_[1] = Vec3f(-1.f, -2.f, -3.f);
        data3f_[2] = Vec3f(2.f, 1.f, -3.f); data3f_[3] = Vec3f(-2.f, -1.f, 3.f);

        // 10 points in 2D forming a "strip" :::::
        data2d_[0] = Vec2d(-2., 1.); data2d_[1] = Vec2d(-1., 1.);
        data2d_[2] = Vec2d(0., 1.);  data2d_[3] = Vec2d(1., 1.);
        data2d_[4] = Vec2d(2., 1.);
        data2d_[5] = Vec2d(-2., 2.); data2d_[6] = Vec2d(-1., 2.);
        data2d_[7] = Vec2d(0., 2.);  data2d_[8] = Vec2d(1., 2.);
        data2d_[9] = Vec2d(2., 2.);

        data1d_[0] = 1.2; data1d_[1] = 2.3; data1d_[2] = 3.4; data1d_[3] = 4.5;
    }

    // Data samples for tests.
    Vecs3f data3f_;
    Vecs2d data2d_;
    Vecs2i data2i_;
    Vecs3f empty1d_;
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

    // An exception should be throw by mean() iff the set is empty.
    EXPECT_THROW(mean(empty1d_), ex_t);
}

TEST_F(MathTest, PCA)
{
    typedef blas::PCA<double, 2> PCAEngine;
    PCAEngine pca;
    PCAEngine::Result result = pca(data2d_);

    PCAEngine::EigenValues eigenvalues = result.get<0>();
    EXPECT_DOUBLE_EQ(0.25, eigenvalues[0]);
    EXPECT_DOUBLE_EQ(2., eigenvalues[1]);

    PCAEngine::EigenVector eigenvector1 = result.get<1>()[0];
    EXPECT_DOUBLE_EQ(0., eigenvector1[0]);
    EXPECT_DOUBLE_EQ(1., eigenvector1[1]);

    PCAEngine::EigenVector eigenvector2 = result.get<1>()[1];
    EXPECT_DOUBLE_EQ(1., eigenvector2[0]);
    EXPECT_DOUBLE_EQ(0., eigenvector2[1]);
}
