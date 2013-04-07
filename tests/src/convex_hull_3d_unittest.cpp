#include <gtest/gtest.h>

#include "bo/methods/convex_hull_3d.hpp"

using namespace bo::methods::surfaces;

// "Test fixture" using base class form GTEST.
class ConvexHull3DTest: public testing::Test
{
protected:

    ConvexHull3DTest()
    { }

    virtual void SetUp()
    {
    }
};

TEST_F(ConvexHull3DTest, UnitCubeHull)
{
    typedef IncrementalConvexHull3D<float> Hull;

    Hull::Points3D points;
    points.push_back(Hull::Point3D(0, 0, 0));
    points.push_back(Hull::Point3D(1, 0, 0));
    points.push_back(Hull::Point3D(0, 1, 0));
    points.push_back(Hull::Point3D(0, 0, 1));
    points.push_back(Hull::Point3D(1, 1, 0));
    points.push_back(Hull::Point3D(0, 1, 1));
    points.push_back(Hull::Point3D(1, 0, 1));
    points.push_back(Hull::Point3D(1, 1, 1));

    Hull conhul3d(points);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 12U);

    float volume = conhul3d.get_volume();

    EXPECT_EQ(volume, 1.0f);
}
