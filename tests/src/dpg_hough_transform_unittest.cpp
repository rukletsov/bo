
#include <cstdlib>
#include <ctime>
#include <gtest/gtest.h>

#include "bo/recognition/dpg_hough_transform_2d.hpp"

using namespace bo;
using namespace bo::recognition;
using namespace bo::recognition::detail;

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

        grid_size_ = Space<float>::Point4D(0.1f, 0.1f, 0.1f, 0.1f);


        // Space-line intersection data.
        pb1_ = Space<float>::Point4D(0.f, 0.f, 0.f, 0.f);
        pb2_ = Space<float>::Point4D(1.f, 1.f, 1.f, 1.f);
        pb3_ = Space<float>::Point4D(1.f, 0.f, 0.f, 0.f);
        pb4_ = Space<float>::Point4D(1.f, 0.f, 1.f, 0.f);
        pb5_ = Space<float>::Point4D(-3.f, -3.f, -3.f, -3.f);
        pb6_ = Space<float>::Point4D(0.5f, 0.5f, 0.5f, 0.5f);
        pb7_ = Space<float>::Point4D(0.1f, 0.f, 0.f, 0.f);
        pb8_ = Space<float>::Point4D(0.4f, 1.f, 0.5f, 0.f);
        space1_ = Space<float>(Space<float>::Box4D(pb1_, pb2_),
                               Space<float>::Size4D(10, 10, 10, 10),
                               1, 1, 0);

        space1x_ = Space<float>(Space<float>::Box4D(Space<float>::Point4D(-1.f, -1.f, -1.f, -1.f),
                                                    Space<float>::Point4D(1.f, 1.f, 1.f, 1.f)));
        pb8a_ = Space<float>::Point4D(1.f, 0.f, 1.f, 0.f);
        pb8b_ = Space<float>::Point4D(0.86f, 0.51f, 1.57f, -0.47f);

        // Space subdivision.
        pb9_ = Space<float>::Point4D(10.f, 10.f, 10.f, 10.f);
        space2_ = Space<float>(Space<float>::Box4D(pb1_, pb9_),
                               Space<float>::Size4D(10, 10, 10, 10),
                               2, 2, 0);


        // Simple object detection.
        // Create a spirale.
        float pi = boost::math::constants::pi<float>();
        for (int i = 0; i < 360; i += 20)
        {
            float angle = i * pi / 180;

            float x = std::cos(angle);
            float y = std::sin(angle);

            DualPointGHT<float>::Feature f(DualPointGHT<float>::Point2D(angle * x, angle * y), DualPointGHT<float>::Point2D(-y, x));
            model_.push_back(f);
        }

        // Create reference points.
        ref3_ = DualPointGHT<float>::Point2D(-1.0f, 0.f);
        ref4_ = DualPointGHT<float>::Point2D(1.0f, 0.f);
        // Object recognition area.
        bbox_ = DualPointGHT<float>::SearchArea(DualPointGHT<float>::Point2D(-2.f, -2.f),
                                                    DualPointGHT<float>::Point2D(2.f, 2.f));


        pb10_ = Space<float>::Point4D(1.f, 1.f, 1.f, 0.9f);



        space3_1_ = Space<float>(Space<float>::Box4D(pb1_, pb2_),
                                 Space<float>::Size4D(3, 3, 3, 3),
                                 1, 1, 0);
        space3_2_ = space3_1_;
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
    Space<float>::Point4D pb8a_;
    Space<float>::Point4D pb8b_;
    Space<float> space1_;
    Space<float> space1x_;
    Space<float> space2_;
    Space<float>::Point4D grid_size_;
    Space<float>::Point4D pb10_;

    DualPointGHT<float>::Features model_;
    DualPointGHT<float>::SearchArea bbox_;
    DualPointGHT<float>::Point2D ref3_;
    DualPointGHT<float>::Point2D ref4_;

    Space<float> space3_1_;
    Space<float> space3_2_;
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
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb1_, pb2_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 2.f);

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb1_, pb3_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 1.f);

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb1_, pb4_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), std::sqrt(2.f));

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb5_, pb2_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 2.f);

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb5_, pb3_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 0.f);

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb6_, pb7_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), 1.0f);

    space1_.reset_votes();
    space1_.vote(space1_.intersect(Space<float>::Line4D(pb7_, pb8_)));

    EXPECT_FLOAT_EQ(space1_.get_votes(), std::sqrt(1 + 0.4f * 0.4f + 0.5f * 0.5f));
}

TEST_F(HoughTransformTest, ComplexSpaceLineIntersection)
{
    space1x_.reset_votes();
    space1x_.vote(space1x_.intersect(Space<float>::Line4D(pb8a_, pb8b_)));

    EXPECT_NEAR(space1x_.get_votes(), 2.44, 0.006);


    // Projection in one dimension is finite null, in another dimension is
    // infinite null.
    space1_.reset_votes();
    Space<float>::Line4D ln(pb5_, pb3_);
    Space<float>::Segment4D sg(ln, DualPointGHT<float>::Point2D(4.f, 5.f));
    space1_.vote_unit(space1_.intersect(sg));

    EXPECT_EQ(space1_.get_votes(), 0);
}

TEST_F(HoughTransformTest, DescreteSpaceLineIntersection)
{
    // Parallel.
    space1_.reset_votes();
    Space<float>::Line4D ln(pb1_, pb3_);
    Space<float>::Segment4D sg(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space1_.vote_descrete(sg);

    EXPECT_EQ(space1_.get_votes(), 10);

    // Diagonal.
    space1_.reset_votes();
    ln = Space<float>::Line4D(pb5_, pb2_);
    sg = Space<float>::Segment4D(ln, DualPointGHT<float>::Point2D(2.f, 4.f));
    space1_.vote_descrete(sg);

    EXPECT_EQ(space1_.get_votes(), 11);

    // Irregular.
    space1_.reset_votes();
    ln = Space<float>::Line4D(pb1_, pb10_);
    sg = Space<float>::Segment4D(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space1_.vote_descrete(sg);

    EXPECT_EQ(space1_.get_votes(), 18);
}


TEST_F(HoughTransformTest, TaxicabSpaceLineIntersection)
{
    // Parallel.
    space1_.reset_votes();
    Space<float>::Line4D ln(pb1_, pb3_);
    Space<float>::Segment4D sg(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space1_.vote_taxicab(sg);

    EXPECT_EQ(space1_.get_votes(), 10);
}

TEST_F(HoughTransformTest, MaxnormSpaceLineIntersection)
{
    // Parallel.
    space1_.reset_votes();
    Space<float>::Line4D ln(pb1_, pb3_);
    Space<float>::Segment4D sg(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space1_.vote_maxnorm(sg);

    EXPECT_EQ(space1_.get_votes(), 10);

    // Diagonal.
    space1_.reset_votes();
    ln = Space<float>::Line4D(pb5_, pb2_);
    sg = Space<float>::Segment4D(ln, DualPointGHT<float>::Point2D(2.f, 4.f));
    space1_.vote_maxnorm(sg);

    EXPECT_EQ(space1_.get_votes(), 10);

    // Irregular.
    space1_.reset_votes();
    ln = Space<float>::Line4D(pb1_, pb10_);
    sg = Space<float>::Segment4D(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space1_.vote_maxnorm(sg);

    EXPECT_EQ(space1_.get_votes(), 10);
}

TEST_F(HoughTransformTest, SpacePlaneIntersection)
{
    Hyperplane4D<float> pl(pb6_, pb3_);
    Space<float>::Points4D p = space1_.intersect(pl);
    EXPECT_EQ(p.size(), 8U);
}

Hyperplane4D<float>::Point4D random_vec()
{
    Hyperplane4D<float>::Point4D p;

    for (int i = 0; i < 4; ++i)
    {
        p[i] = 100 * (std::rand() + 0.0f) / RAND_MAX;
    }

    return p;
}

TEST_F(HoughTransformTest, PlaneIntersectionDecomposition)
{
    for (std::size_t i = 0; i < 100; ++i)
    {
        Hyperplane4D<float> plane(random_vec(), random_vec());

        Hyperplane4D<float>::Point4D pd = plane.decompose(plane.point());

        EXPECT_NEAR(pd[3], 0, 0.001);

        Hyperplane4D<float>::Point4D q1 = random_vec();
        Hyperplane4D<float>::Point4D q2 = random_vec();
        float t;
        if (plane.intersect(q1, q2, t))
        {
            Hyperplane4D<float>::Point4D inters = q1 + t * (q2 - q1);

            pd = plane.decompose(inters);

            EXPECT_NEAR(pd[3], 0, 0.001);
        }

    }
}

TEST_F(HoughTransformTest, SpaceVolumeSubdivision)
{
    Space<float> s = Space<float>(Space<float>::Box4D(pb1_, pb2_),
                           Space<float>::Size4D(2, 2, 2, 2),
                           1, 1, 0);

    Hyperplane4D<float> pl(pb6_ / 2, pb3_);
    s.vote(pl);

    EXPECT_NEAR(s.get_votes(), 1, 0.001);

    s.subdivide();
    float votes_acc = 0;

    for (Space<float>::Spaces::iterator it = s.get_subspaces().begin();
         it != s.get_subspaces().end(); ++it)
    {
        it->vote(pl);
        votes_acc += it->get_votes();
    }

    EXPECT_NEAR(s.get_votes(), votes_acc, 0.001);
}

TEST_F(HoughTransformTest, SpaceSubdivision)
{
    space2_.subdivide();
    Space<float>::Spaces subs = space2_.get_subspaces();

    EXPECT_EQ(subs.size(), 10000u);
}

TEST_F(HoughTransformTest, SelfDetection)
{
//    // Create dual-point reference.
//    DualPointGHT<float>::Reference ref(ref3_, ref4_);

//    // Create a transformation and encode the feature.
//    DualPointGHT<float> ght(model_, ref, 0.01f);

//    DualPointGHT<float>::ReferenceVotes rv = ght.cross_level_detect(
//                                                model_, 0.7f, Space<float>::Size4D(3, 3, 3, 3),
//                                                6, 0, bbox_, bbox_);

//    DualPointGHT<float>::Reference r = rv.front().first;

//    EXPECT_NEAR(r.first[0], -1, 0.01);
//    EXPECT_NEAR(r.first[1],  0, 0.01);
//    EXPECT_NEAR(r.second[0],  1, 0.01);
//    EXPECT_NEAR(r.second[1],  0, 0.01);
}

TEST_F(HoughTransformTest, ProbabilisticModel)
{
    std::size_t k = 79;
    std::size_t n = 37;

    // Test of the distribution function.
    float f = SubdivisionPolicy<float>::F(k, n, k);

    EXPECT_NEAR(f, 1.0, 0.0001);

    Space<float>::Line4D ln(pb1_, pb3_);
    Space<float>::Segment4D sg(ln, DualPointGHT<float>::Point2D(0.f, 1.f));
    space3_1_.vote_maxnorm(sg);
    space3_1_.vote_maxnorm(sg);
    space3_2_.vote_maxnorm(sg);
    Space<float>::Spaces subc;
    subc.push_back(space3_1_);
    subc.push_back(space3_2_);

    float p = SubdivisionPolicy<float>::p_max_value_in_subcollection(subc, subc);

    // P(max1 > max2) = P(max2 > max1), and
    // P(max1 > max2) + P(max2 > max1) + P(max1 = max2) = 1.
    EXPECT_EQ(2 * p < 1, true);
}
