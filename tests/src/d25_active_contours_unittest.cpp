#include <gtest/gtest.h>

#include "bo/methods/d25_active_contours.hpp"

using namespace bo;
using namespace bo::methods::surfaces;
using namespace bo::methods::surfaces::detail;

// "Test fixture" using base class form GTEST.
class D25ActiveContoursTest: public testing::Test
{
protected:

    D25ActiveContoursTest()
    { }

    virtual void SetUp()
    {
        tr1_1_ = Triangle<Vertex>(Vertex(0, 0, 0), Vertex(1, 0, 0), Vertex(0, 1, 0));
        tr1_2_ = Triangle<Vertex>(Vertex(1, 1, 0), Vertex(1, 0, 0), Vertex(0, 1, 0));
        height1_ = 1.0f;
        angle1_ = 0.5f;

        tr2_ = Triangle<Vertex>(Vertex(0, 0, 1), Vertex(1, 0, 1), Vertex(0, 1, 1));
        height2_1_ = 0.5f;
        height2_2_ = 0.51f;
        angle2_1_ = 0.505458f;
        angle2_2_ = angle2_1_ - 0.01f;
    }

    float height1_;
    float angle1_;
    Triangle<Vertex> tr1_1_;
    Triangle<Vertex> tr1_2_;

    float height2_1_;
    float height2_2_;
    float angle2_1_;
    float angle2_2_;
    Triangle<Vertex> tr2_;
};

TEST_F(D25ActiveContoursTest, TriangularDipyramidIntersection)
{
    TriangularDipyramid tdp1 = TriangularDipyramid::from_triangle_and_height(tr1_1_, height1_);
    TriangularDipyramid tdp2 = TriangularDipyramid::from_triangle_and_height(tr1_2_, height1_);
    EXPECT_EQ(tdp1.intersects(tdp2), false);
    EXPECT_EQ(tdp1.intersects(tdp1), true);

    tdp1 = TriangularDipyramid::from_triangle_and_angle(tr1_1_, angle1_);
    tdp2 = TriangularDipyramid::from_triangle_and_angle(tr1_2_, angle1_);
    EXPECT_EQ(tdp1.intersects(tdp2), false);
    EXPECT_EQ(tdp1.intersects(tdp1), true);

    tdp1 = TriangularDipyramid::from_triangle_and_height(tr1_1_, height2_1_);
    tdp2 = TriangularDipyramid::from_triangle_and_height(tr2_, height2_1_);
    EXPECT_EQ(tdp1.intersects(tdp2), false);

    tdp1 = TriangularDipyramid::from_triangle_and_height(tr1_1_, height2_2_);
    tdp2 = TriangularDipyramid::from_triangle_and_height(tr2_, height2_2_);
    EXPECT_EQ(tdp1.intersects(tdp2), true);

    tdp1 = TriangularDipyramid::from_triangle_and_angle(tr1_1_, angle2_1_);
    tdp2 = TriangularDipyramid::from_triangle_and_angle(tr2_, angle2_1_);
    EXPECT_EQ(tdp1.intersects(tdp2), false);

    tdp1 = TriangularDipyramid::from_triangle_and_angle(tr1_1_, angle2_2_);
    tdp2 = TriangularDipyramid::from_triangle_and_angle(tr2_, angle2_2_);
    EXPECT_EQ(tdp1.intersects(tdp2), true);
}
