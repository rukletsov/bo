
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

    RawImage2DTest(): im1_(400, 300), im2_(5, 1), im3_(1, 7), im_invalid1_(0, 6)
    { }

    virtual void SetUp()
    {
        // Nullify im3_ image.
        memset(im3_.data(), 0, im3_.size() * sizeof(boost::uint8_t));

        // Initialize im2_ image.
        for (std::size_t i = 0; i < im2_.width(); ++i)
            im2_.at(i, 0) = double(i);
    }

    // Variables for all RawImage2D<> tests to use.
    RawImage2D<float> im1_;
    RawImage2D<double> im2_;
    RawImage2D<boost::uint8_t> im3_;
    RawImage2D<float> im_invalid1_;
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
    EXPECT_TRUE(im_invalid1_.is_null());
    EXPECT_EQ(std::size_t(0), im_invalid1_.size());
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

TEST_F(RawImage2DTest, SettersGetters)
{
    // These checks are done implicitly in almost all other tests (as well during
    // SetUp() call). However several tests will be provided here in order to ensure
    // the coverage of the test.
    EXPECT_DOUBLE_EQ(0., im2_(0, 0));
    EXPECT_DOUBLE_EQ(1., im2_.at(1, 0));
    EXPECT_DOUBLE_EQ(3., im2_(3, 0));
    EXPECT_DOUBLE_EQ(4., im2_.at(4, 0));

    im3_(0, 1) = boost::uint8_t(5);
    im3_.at(0, 6) = boost::uint8_t(4);

    std::size_t sum(0);
    for (std::size_t i = 0; i < im3_.size(); ++i)
        sum += im3_.data()[i];

    EXPECT_EQ(std::size_t(9), sum);
}

TEST_F(RawImage2DTest, BoundaryChecks)
{
    typedef std::out_of_range ex_t;

    // An exception should be throw by at() iff the given index is out of range.
    EXPECT_NO_THROW(im2_.at(0, 0));
    EXPECT_NO_THROW(im2_.at(4, 0) = 1.);

    EXPECT_THROW(im2_.at(-1, 0), ex_t);
    EXPECT_THROW(im2_.at(0, -1), ex_t);
    EXPECT_THROW(im2_.at(5, 0) = 1., ex_t);
    EXPECT_THROW(im2_.at(3, 1) = -1., ex_t);

    // For invalid images (image data is NULL), every index is out of range.
    EXPECT_THROW(im_invalid1_.at(0, 0) = 1.f, ex_t);
    EXPECT_THROW(im_invalid1_.at(0, 3), ex_t);
    EXPECT_THROW(im_invalid1_.at(1, 6) = -1.f, ex_t);

    // Operator () doesn't throw exceptions when wrong index is used. Instead a
    // debug assertion is called. Since assertion leads to a special state of a program
    // (usually termination and via abort()), corresponding tests are placed in a
    // special test case for so-called "death tests". However, valid indices should
    // be always processed without any crashes.
    EXPECT_NO_THROW(im3_(0, 0) = boost::uint8_t(1));
    EXPECT_NO_THROW(im3_(0, 2) = boost::uint8_t(0));
    EXPECT_NO_THROW(im3_(0, 6));
}

TEST_F(RawImage2DTest, NullCheck)
{
    // Method is_null() should return true iff the image is invalid. If is_null()
    // returns true, at() must throw an exception for any index.
    EXPECT_FALSE(im1_.is_null());
    EXPECT_FALSE(im2_.is_null());
    EXPECT_FALSE(im3_.is_null());

    EXPECT_TRUE(im_invalid1_.is_null());
    EXPECT_TRUE((RawImage2D<int>()).is_null());
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

    EXPECT_DEATH(im_invalid1_.offset(0, 0), ".*");
    EXPECT_DEATH(im_invalid1_.offset(-1, 1), ".*");
}

TEST_F(RawImage2DDeathTest, BoundaryChecksAssertions)
{
    // Accessing image pixel data by a wrong index should lead to a debug assertion.
    // It is expected to terminate the program in debug mode. Note that at() function
    // should not use assertions (but exceptions) and therefore should not be checked.
    EXPECT_DEATH(im3_(0, -1), ".*");
    EXPECT_DEATH(im3_(-1, 0) = boost::uint8_t(0), ".*");
    EXPECT_DEATH(im3_(0, 7) = boost::uint8_t(1), ".*");
    EXPECT_DEATH(im3_(1, 2), ".*");

    // For invalid images (image data is NULL), every index is out of range.
    EXPECT_DEATH(im_invalid1_(0, -1), ".*");
    EXPECT_DEATH(im_invalid1_(im_invalid1_.width() - 1, im_invalid1_.height() - 1)
                 = 0.f, ".*");
    EXPECT_DEATH(im_invalid1_(-1, 6), ".*");
}

#endif
