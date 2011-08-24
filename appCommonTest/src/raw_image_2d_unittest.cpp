
#include <cstddef>
#include <boost/cstdint.hpp>
#include <gtest/gtest.h>

#include "common/raw_image_2d.hpp"
#include "debug_alloc.hpp"

using namespace common;


// Create a so-called "text fixture" using base class form GTEST.
class RawImage2DTest: public testing::Test
{
protected:

    RawImage2DTest(): im1_(400, 300), im2_(5, 1), im3_(1, 7), im4_(0, 6)
    { }

    // Variables for all RawImage2D<> tests to use.
    RawImage2D<float> im1_;
    RawImage2D<double> im2_;
    RawImage2D<boost::uint8_t> im3_;
    RawImage2D<float> im4_;
};

TEST_F(RawImage2DTest, DefaultConstructor)
{
    // Default constructor should allocate no memory for the image.
    RawImage2D<double> image;
    EXPECT_EQ(std::size_t(0), image.width());
    EXPECT_EQ(std::size_t(0), image.height());
    EXPECT_EQ(NULL, image.data());
}

TEST_F(RawImage2DTest, SizeConstructor)
{
    // An uninitialized array shoud be allocated for image data. Its size should
    // correspond with given size during cunstruction.
    EXPECT_EQ(std::size_t(400), im1_.width());
    EXPECT_EQ(std::size_t(300), im1_.height());

    EXPECT_EQ(std::size_t(1), im3_.width());
    EXPECT_EQ(std::size_t(1), im2_.height());

    // TODO: check if the last pixel (x + width * y - 1) belongs to the image and
    // the next byte doesn't.

    // An image with one of the dimensions equals zero should not be created.
    EXPECT_TRUE(im4_.is_null());
    EXPECT_EQ(std::size_t(0), im4_.size());
}

TEST_F(RawImage2DTest, OffsetCalculation)
{
    // Method offset() doesn't throw exceptions when wrong index is used. Instead a
    // debug assertion is used. Since assertion leads to a special state of a program
    // (usually termination and via abort()), corresponding tests are placed in a
    // special test case for so-called "death tests". However, valid indices should
    // be always processed without any crashes.
    EXPECT_EQ(std::size_t(0), im1_.offset(0, 0));
    EXPECT_EQ(std::size_t(im1_.width()), im1_.offset(0, 1));
    EXPECT_EQ(std::size_t(12824), im1_.offset(24, 32));
    EXPECT_EQ(std::size_t(im1_.size() - 1), im1_.offset(399, 299));

    EXPECT_EQ(std::size_t(4), im2_.offset(4, 0));
}


// An aliased fixture for so-called "death tests".
typedef RawImage2DTest RawImage2DDeathTest;

#ifdef _DEBUG

TEST_F(RawImage2DDeathTest, OffsetAssertions)
{
    // Calculating offset for invalid index should lead to a debug assertion, which
    // expected to entail a program termination in debug mode (at least for console
    // applications).
    EXPECT_DEATH(im1_.offset(400, 299), ".*");
    EXPECT_DEATH(im1_.offset(399, 300), ".*");
    EXPECT_DEATH(im1_.offset(-1, 0), ".*");
    EXPECT_DEATH(im1_.offset(399, -1), ".*");

    EXPECT_DEATH(im2_.offset(4, 1), ".*");
    EXPECT_DEATH(im3_.offset(1, 1), ".*");

    EXPECT_DEATH(im4_.offset(0, 0), ".*");
    EXPECT_DEATH(im4_.offset(-1, 1), ".*");
}

#endif
