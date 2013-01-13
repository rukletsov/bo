#include <gtest/gtest.h>

#include "bo/methods/dpg_hough_transform_2d.hpp"

using namespace bo;
using namespace bo::methods::recognition;


// "Test fixture" using base class form GTEST.
class HoughTransformTest: public testing::Test
{
protected:

    HoughTransformTest(): ref1_(0.f, 0.f), ref2_(1.f, 1.f)
    { }

    virtual void SetUp()
    {
        c1_ = DualPointGHT<float>::Point2D(1.f,0.f);
        tang1_ = DualPointGHT<float>::Point2D(1.f, 1.f);

        c2_ = DualPointGHT<float>::Point2D(5.f, 5.f);
        tang2_ = DualPointGHT<float>::Point2D(-1.f, -1.f);
    }

    DualPointGHT<float>::Point2D ref1_;
    DualPointGHT<float>::Point2D ref2_;

    DualPointGHT<float>::Point2D c1_;
    DualPointGHT<float>::Point2D tang1_;
    DualPointGHT<float>::Point2D c2_;
    DualPointGHT<float>::Point2D tang2_;
};

TEST_F(HoughTransformTest, EncodeReconstruct)
{
    // Create a collection of model features.
    DualPointGHT<float>::Feature f1(c1_, tang1_);
    DualPointGHT<float>::Feature f2(c2_, tang2_);
    DualPointGHT<float>::Features fs;
    fs.push_back(f1);
    fs.push_back(f2);

    // Create dual-point reference.
    DualPointGHT<float>::Reference ref(ref1_, ref2_);

    // Create a transformation and encode the feature.
    DualPointGHT<float> ght(fs, ref);

    // Decode the feature back.
    DualPointGHT<float>::Points2D ps = ght.reconstruct(ref);

    EXPECT_EQ(ps.size(), 2U);

    EXPECT_NEAR(ps.front()[0], f1.first[0], 0.001f);
    EXPECT_NEAR(ps.front()[1], f1.first[1], 0.001f);

    EXPECT_NEAR(ps.back()[0], f2.first[0], 0.001f);
    EXPECT_NEAR(ps.back()[1], f2.first[1], 0.001f);

}
