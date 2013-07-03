
#include <gtest/gtest.h>

#include "bo/math/topology.hpp"

using namespace bo;
using namespace bo::math;

TEST(TopologyTest, OrthotopeEdgeCount)
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
