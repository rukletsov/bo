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
    EXPECT_EQ((OrthotopeTopology<float, 1>::edges().size()),
              (OrthotopeTopology<float, 1>::edge_count()));

    EXPECT_EQ((OrthotopeTopology<float, 2>::edges().size()),
              (OrthotopeTopology<float, 2>::edge_count()));

    EXPECT_EQ((OrthotopeTopology<float, 3>::edges().size()),
              (OrthotopeTopology<float, 3>::edge_count()));

    EXPECT_EQ((OrthotopeTopology<float, 4>::edges().size()),
              (OrthotopeTopology<float, 4>::edge_count()));

    EXPECT_EQ((OrthotopeTopology<float, 5>::edges().size()),
              (OrthotopeTopology<float, 5>::edge_count()));
}
