
#include <gtest/gtest.h>

#include "bo/distances/distances_3d.hpp"
#include "debug_alloc.hpp"

using namespace bo;
using namespace bo::distances;

// Create a so-called "text fixture".
class GeometryTest: public testing::Test
{
public:
    typedef Vector<float, 3> Vec3f;

protected:

    GeometryTest(): plane1_origin_(0, 0, 0), plane1_norm_(0, 0, 1)
    { }

    Vec3f plane1_origin_;
    Vec3f plane1_norm_;
};


TEST_F(GeometryTest, PlaneProjection)
{
    Vec3f proj1 = project_point_onto_plane(Vec3f(1.f, 1.f, 1.f), plane1_origin_, plane1_norm_);
    EXPECT_FLOAT_EQ(1.f, proj1[0]);
    EXPECT_FLOAT_EQ(1.f, proj1[1]);
    EXPECT_FLOAT_EQ(0.f, proj1[2]);
}
