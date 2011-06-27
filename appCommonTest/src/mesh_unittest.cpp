
#include <gtest/gtest.h>

#include "common/mesh.hpp"

using namespace common;


// Create a so-called "text fixture" using base class form GTEST.
class MeshTest: public testing::Test
{
protected:

    MeshTest(): mesh1_(10)
    { }

    virtual void SetUp()
    { }

    // Create a test mesh object with 10 vertices.
    Mesh mesh1_;
};


TEST_F(MeshTest, SimpleTest)
{
    EXPECT_EQ(0, mesh1_.get_all_faces().size());
}
