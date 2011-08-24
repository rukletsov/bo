
#include <cstddef>
#include <gtest/gtest.h>

#include "common/raw_image_2d.hpp"

using namespace common;


// Create a so-called "text fixture" using base class form GTEST.
class RawImage2DTest: public testing::Test
{
protected:

    RawImage2DTest(): im1_(400, 300)
    { }

    // Variables for all Vector<> tests to use.
    RawImage2D<float> im1_;
};

TEST_F(RawImage2DTest, DefaultConstructor)
{
    // Default constructor should allocate no memory for the image.
    RawImage2D<double> image;
    EXPECT_EQ(std::size_t(0), image.width());
    EXPECT_EQ(std::size_t(0), image.height());
    EXPECT_EQ(NULL, image.data());
}
