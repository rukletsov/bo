#include <gtest/gtest.h>

#include "bo/methods/dpg_hough_transform_2d.hpp"

using namespace bo;
using namespace bo::methods::recognition;
using namespace bo::methods::recognition::detail;



// "Test fixture" using base class form GTEST.
class HoughTransformTest: public testing::Test
{
protected:

    HoughTransformTest(): ref1_(0.f, 0.f), ref2_(1.f, 1.f)
    { }

    virtual void SetUp()
    {
        // Encode-reconstruct data.
        c1_ = DualPointGHT<float>::Point2D(1.f,0.f);
        tang1_ = DualPointGHT<float>::Point2D(1.f, 1.f);

        c2_ = DualPointGHT<float>::Point2D(5.f, 5.f);
        tang2_ = DualPointGHT<float>::Point2D(-1.f, -1.f);


        // Space-line intersection data.
        pb1_ = Space<float>::Point4D(0.f, 0.f, 0.f, 0.f);
        pb2_ = Space<float>::Point4D(1.f, 1.f, 1.f, 1.f);
        pb3_ = Space<float>::Point4D(1.f, 0.f, 0.f, 0.f);
        pb4_ = Space<float>::Point4D(1.f, 0.f, 1.f, 0.f);
        pb5_ = Space<float>::Point4D(-3.f, -3.f, -3.f, -3.f);
        pb6_ = Space<float>::Point4D(0.5f, 0.5f, 0.5f, 0.5f);
        pb7_ = Space<float>::Point4D(0.1f, 0.f, 0.f, 0.f);
        pb8_ = Space<float>::Point4D(0.4f, 1.f, 0.5f, 0.f);
        space1_ = Space<float>(Space<float>::Box4D(pb1_, pb2_));

        // Space subdivision.
        pb9_ = Space<float>::Point4D(10.f, 10.f, 10.f, 10.f);
        space2_ = Space<float>(Space<float>::Box4D(pb1_, pb9_));

        // Object detection.
        // Generate a circle model with radius 1 and center (0,0).
        for (float a = 0; a < boost::math::constants::pi<float>(); a += 0.1)
        {
            DualPointGHT<float>::Point2D p(std::cos(a), std::sin(a));
            DualPointGHT<float>::Point2D tan(-std::sin(a), std::cos(a));

            model.push_back( DualPointGHT<float>::Feature(p, tan));
        }
        // Generate a circle object with radius 1 and center (2, 2).
        for (float a = 0; a < boost::math::constants::pi<float>(); a += 0.1)
        {
            DualPointGHT<float>::Point2D p(2 + std::cos(a), 2 + std::sin(a));
            DualPointGHT<float>::Point2D tan(-std::sin(a), std::cos(a));

            object.push_back( DualPointGHT<float>::Feature(p, tan));
        }
        // Create reference points.
        ref3_ = DualPointGHT<float>::Point2D(-0.5f,0.f);
        ref4_ = DualPointGHT<float>::Point2D(0.5f,0.f);
        // Object recognition area.
        bbox = DualPointGHT<float>::RecognitionArea(DualPointGHT<float>::Point2D(-1.f, -1.f),
                                                    DualPointGHT<float>::Point2D(3.f, 3.f));
    }

    DualPointGHT<float>::Point2D ref1_;
    DualPointGHT<float>::Point2D ref2_;

    DualPointGHT<float>::Point2D c1_;
    DualPointGHT<float>::Point2D tang1_;
    DualPointGHT<float>::Point2D c2_;
    DualPointGHT<float>::Point2D tang2_;


    Space<float>::Point4D pb1_;
    Space<float>::Point4D pb2_;
    Space<float>::Point4D pb3_;
    Space<float>::Point4D pb4_;
    Space<float>::Point4D pb5_;
    Space<float>::Point4D pb6_;
    Space<float>::Point4D pb7_;
    Space<float>::Point4D pb8_;
    Space<float>::Point4D pb9_;
    Space<float> space1_;
    Space<float> space2_;

    DualPointGHT<float>::Features model;
    DualPointGHT<float>::Features object;
    DualPointGHT<float>::RecognitionArea bbox;
    DualPointGHT<float>::Point2D ref3_;
    DualPointGHT<float>::Point2D ref4_;
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

TEST_F(HoughTransformTest, SpaceLineIntersection)
{
    space1_.intersect(Space<float>::Line4D(pb1_, pb2_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 2.f);

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb1_, pb3_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 1.f);

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb1_, pb4_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), std::sqrt(2.f));

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb5_, pb2_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 2.f);

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb5_, pb3_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 0.f);

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb6_, pb7_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 0.5f);

    space1_.reset_votes();
    space1_.intersect(Space<float>::Line4D(pb7_, pb8_));

    EXPECT_FLOAT_EQ(space1_.get_votes(), std::sqrt(1 + 0.4f * 0.4f + 0.5f * 0.5f));

}

TEST_F(HoughTransformTest, SpaceSubdivision)
{
    std::size_t n = 5;

    space2_.subdivide(n);
    detail::Space<float>::Spaces subs = space2_.get_subspaces();

    EXPECT_EQ(subs.size(), n * n * n * n);
}

TEST_F(HoughTransformTest, SelfDetection)
{
    // Create dual-point reference.
    DualPointGHT<float>::Reference ref(ref3_, ref4_);

    // Create a transformation and encode the feature.
    DualPointGHT<float> ght(model, ref);

    DualPointGHT<float>::ReferenceVotes rv = ght.fast_detect(model, 0.9f, 2, 1, bbox);

    DualPointGHT<float>::Reference r = rv.front().first;

    EXPECT_EQ(r.first.x() > 0 && r.first.x() < 1 &&
              r.first.y() > 0 && r.first.y() < 1 &&
              r.second.x() > 0 && r.second.x() < 1 &&
              r.second.y() > 0 && r.second.y() < 1, true);

    int iii = 0;
}
