#include <gtest/gtest.h>

#include "bo/topology.hpp"

using namespace bo;
using namespace bo::topology;

// "Test fixture" using base class form GTEST.
class TopologyTest: public testing::Test
{
protected:

    TopologyTest()
    { }

    virtual void SetUp()
    {
    }
};

TEST_F(TopologyTest, OrthotopeEdgeCount)
{
    EXPECT_EQ((OrthotopeGeometry<float, 1>::edges().size()),
              (OrthotopeGeometry<float, 1>::edge_count()));

    EXPECT_EQ((OrthotopeGeometry<float, 2>::edges().size()),
              (OrthotopeGeometry<float, 2>::edge_count()));

    EXPECT_EQ((OrthotopeGeometry<float, 3>::edges().size()),
              (OrthotopeGeometry<float, 3>::edge_count()));

    EXPECT_EQ((OrthotopeGeometry<float, 4>::edges().size()),
              (OrthotopeGeometry<float, 4>::edge_count()));

    EXPECT_EQ((OrthotopeGeometry<float, 5>::edges().size()),
              (OrthotopeGeometry<float, 5>::edge_count()));
}
