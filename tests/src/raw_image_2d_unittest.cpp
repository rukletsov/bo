
#include <cstddef>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <gtest/gtest.h>

#include "bo/core/raw_image_2d.hpp"
#include "debug_alloc.hpp"

using namespace bo;

// Create a so-called "text fixture".
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

    template <typename T>
    static RawImage2D<T>* bw_stripes(std::size_t width, std::size_t height)
    {
        RawImage2D<T>* retvalue = new RawImage2D<T>(width, height);

        // Fill image data with black and white stripes.
        for (std::size_t col = 0; col < retvalue->width(); ++col)
        {
            T val = T(0);
            if ((col % 4 == 2) || (col % 4 == 3))
                val = T(1);
            for (std::size_t row = 0; row < retvalue->height(); ++row)
                retvalue->operator ()(col, row) = val;
        }

        return retvalue;
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

TEST_F(RawImage2DTest, CopyConstructor)
{
    // Copy c-tor should create a wrapper object for the same image data.
    RawImage2D<float> image(im1_);

    // Dimensions should be the same.
    EXPECT_EQ(im1_.size(), image.size());
    EXPECT_EQ(im1_.width(), image.width());
    EXPECT_EQ(im1_.height(), image.height());

    // Changing data through one object should be visible through others.
    image.at(1, 1) = 2.f;
    EXPECT_FLOAT_EQ(2.f, im1_.at(1, 1));

    im1_(14, 78) *= 2;
    EXPECT_FLOAT_EQ(im1_(14, 78), image(14, 78));

    // One more copy.
    RawImage2D<float> image2 = image;

    // Image data should be the same.
    EXPECT_EQ(im1_.data(), image2.data());
    int diff = memcmp(image2.data(), image.data(), im1_.size());
    EXPECT_EQ(int(0), diff);
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

TEST_F(RawImage2DTest, Cloning)
{
    // Method clone() should create a deep copy of the object.
    RawImage2D<float> image = im1_.clone();

    // Dimensions should be the same.
    EXPECT_EQ(im1_.size(), image.size());
    EXPECT_EQ(im1_.width(), image.width());
    EXPECT_EQ(im1_.height(), image.height());

    // Right after the cloning image data should be the same in both objects.
    int diff = memcmp(image.data(), im1_.data(), im1_.size());
    EXPECT_EQ(int(0), diff);

    // But objects should have separate image data.
    RawImage2D<float>::Index index(12, 248);
    im1_(index.first, index.second) += 10.f;
    EXPECT_PRED_FORMAT2(::testing::FloatLE, image(index.first, index.second),
                        im1_(index.first, index.second));
    EXPECT_FLOAT_EQ(im1_(index.first, index.second),
                    image(index.first, index.second) + 10.f);
}

TEST_F(RawImage2DTest, AssignmentOperator)
{
    RawImage2D<boost::uint8_t> image(5, 5);

    // Assignment operator makes a shallow copy of the object.
    image = im3_;

    // Dimensions should be the same.
    EXPECT_EQ(im3_.size(), image.size());
    EXPECT_EQ(im3_.width(), image.width());
    EXPECT_EQ(im3_.height(), image.height());

    // As with copy c-tor, image data is shared among shallow copies.
    image(0, 0) = 74;
    im3_(0, 6) = 9;
    int diff = memcmp(image.data(), im3_.data(), image.size());
    EXPECT_EQ(int(0), diff);
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

    EXPECT_EQ(std::size_t(4), im2_.offset(RawImage2D<double>::Index(4, 0)));
    EXPECT_EQ(std::size_t(0), im3_.offset(RawImage2D<double>::Index(0, 0)));

    EXPECT_EQ(im3_.offset(0, 6), im3_.offset(RawImage2D<boost::uint8_t>::Index(0, 6)));
}

TEST_F(RawImage2DTest, IndexCalculation)
{
    // Method index() doesn't throw exceptions when wrong offset is used. It behaves
    // like offset() if a bad index is used.
    RawImage2D<float>::Index index;

    index = im1_.index(0);
    EXPECT_EQ(std::size_t(0), index.first);
    EXPECT_EQ(std::size_t(0), index.second);

    index = im1_.index(61801); // [201, 154]
    EXPECT_EQ(std::size_t(201), index.first);
    EXPECT_EQ(std::size_t(154), index.second);

    index = im1_.index(im1_.offset(18, 237)); // [18, 237]
    EXPECT_EQ(std::size_t(18), index.first);
    EXPECT_EQ(std::size_t(237), index.second);

    index = im1_.index(23544);
    EXPECT_EQ(std::size_t(23544), im1_.offset(index));

    index = im1_.index(im1_.size() - 1);
    EXPECT_EQ(std::size_t(im1_.width() - 1), index.first);
    EXPECT_EQ(std::size_t(im1_.height() - 1), index.second);
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

    // Check the return value of data() method.
    EXPECT_TRUE(NULL != im1_.data());
    EXPECT_TRUE(NULL != im2_.data());
    EXPECT_TRUE(NULL != im3_.data());

    EXPECT_EQ(NULL, im_invalid1_.data());
    EXPECT_EQ(NULL, (RawImage2D<int>()).data());
}

TEST_F(RawImage2DTest, SizeCheck)
{
    // Method size() should not throw and should always be zero for invalid images.
    // It should correlate with width() and height() methods as well.
    EXPECT_EQ(std::size_t(5), im2_.size());
    EXPECT_EQ(std::size_t(0), im_invalid1_.size());

    EXPECT_EQ(im3_.size(), im3_.width() * im3_.height());
    EXPECT_EQ(im_invalid1_.size(), im_invalid1_.width() * im_invalid1_.height());
}

TEST_F(RawImage2DTest, Fill)
{
    im1_.fill(-1.f);
    EXPECT_EQ(*std::min_element(im1_.data(), im1_.data() + im1_.size()),
              *std::max_element(im1_.data(), im1_.data() + im1_.size()));
    EXPECT_EQ(-1.f, *std::min_element(im1_.data(), im1_.data() + im1_.size()));

    im_invalid1_.fill(1.f);
    EXPECT_EQ(std::size_t(0), im_invalid1_.size());
}

// TODO: provide tests for other RawImage2D functions.


// An aliased fixture for so-called "death tests".
typedef RawImage2DTest RawImage2DDeathTest;

#ifdef _DEBUG

TEST_F(RawImage2DDeathTest, OffsetAssertions)
{
    // Calculating offset for invalid index should lead to a debug assertion, which
    // expected to entail a program termination in debug mode (at least for console
    // applications).
    EXPECT_DEATH(im1_.offset(RawImage2D<float>::Index(400, 299)), ".*");
    EXPECT_DEATH(im1_.offset(RawImage2D<float>::Index(399, 300)), ".*");
    EXPECT_DEATH(im1_.offset(-1, 0), ".*");
    EXPECT_DEATH(im1_.offset(399, -1), ".*");

    EXPECT_DEATH(im2_.offset(4, 1), ".*");
    EXPECT_DEATH(im3_.offset(1, 1), ".*");

    EXPECT_DEATH(im_invalid1_.offset(0, 0), ".*");
    EXPECT_DEATH(im_invalid1_.offset(RawImage2D<float>::Index(-1, 1)), ".*");
}

TEST_F(RawImage2DDeathTest, IndexAssertions)
{
    // Expected to work very similar to offset() method.
    EXPECT_DEATH(im1_.index(-1), ".*");
    EXPECT_DEATH(im3_.index(im3_.size()), ".*");

    EXPECT_DEATH(im_invalid1_.index(0), ".*");
    EXPECT_DEATH(im_invalid1_.index(1), ".*");
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
