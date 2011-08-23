
/******************************************************************************

    raw_image.hpp, v 1.0.3 2011.03.14

    2D image class. OpenCV library can be used for IO. Be advised that
    different versions of OpenCV can have different interface and influence
    RawImage class. For this reason OpenCV 2.0, 2.1 or 2.2 are recommended.

    Copyright (c) 2010, 2011
    Alexander Rukletsov <rukletsov@gmail.com>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
	    notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
	    notice, this list of conditions and the following disclaimer in the
	    documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
    OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGE.

*******************************************************************************/

#ifndef RAW_IMAGE_HPP_A9C93511_7D52_457E_9B7A_5CFA9590A8C9_
#define RAW_IMAGE_HPP_A9C93511_7D52_457E_9B7A_5CFA9590A8C9_

#include <vector>
#include <cmath>

#include <boost/cstdint.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/utility.hpp>

// Enable OpenCV usage only when requested by the library user.
#ifdef USE_OPENCV
    // Depending on the version of OpenCV you link to one of the following 
    // headers should be used. If your OpenCV version is 2.2, define OPENCV_2_2 
    // symbol before including this file. If your OpenCV version is either 2.0 
    // or 2.1, include this file without defining any symbols.
    #ifdef OPENCV_2_2
        #include <opencv2/opencv.hpp>
    #else
        #include <opencv/cv.h>
    #endif
#endif

namespace common {

template <typename ValType>
class RawImage : boost::noncopyable
{

//    typedef std::vector<Pixels> PixelMatrix;
    typedef boost::scoped_ptr<ValType> ImageDataPtr;

public:
    typedef ValType* iterator;
    typedef const ValType* const_iterator;
    typedef ValType& reference;
    typedef const ValType& const_reference;

    typedef std::vector<ValType> Pixels;
    typedef std::pair<std::size_t, std::size_t> Index;
    typedef std::vector<Index> Indices;

public:
    RawImage();
    RawImage(std::size_t width, std::size_t height);

//    RawImage(const PixelMatrix& image);

    // Convertion functions from and to OpenCV format are available on demand.
#ifdef USE_OPENCV
    template <typename T> static
    RawImage<ValType> from_cvmat(const cv::Mat& image);

    cv::Mat to_cvmat() const;
#endif

//    PixelMatrix raw() const;
    ValType* data() const;

    std::size_t size() const;

//    Pixels& operator[](std::size_t row);
//    const Pixels& operator[](std::size_t row) const;

    ValType& at(std::size_t col, std::size_t row);

    Pixels get_neighbour_values(std::size_t row, std::size_t col) const;
    Indices get_neighbours(std::size_t row, std::size_t col) const;

protected:
    double av_dist(std::size_t row, std::size_t col) const;
    double std_devia(std::size_t row, std::size_t col) const;

private:
    std::size_t width_;
    std::size_t height_;

    ImageDataPtr image_;
};


template <typename ValType>
RawImage<T>::RawImage() : width_(0), height_(0), image_(NULL)
{ }

template <typename ValType>
RawImage::RawImage(std::size_t width, std::size_t height) : width_(width),
    height_(height), image_(new ValType[width * height])
{ }

#ifdef USE_OPENCV

// Utility functions for RawImage class with OpenCV usage.
template <typename T> inline 
bool is_bpps_equal_to(const cv::Mat& model)
{
    return 
        (sizeof(T) == model.elemSize());
}

inline
void show_image(const std::string& caption, const cv::Mat& image)
{
    cv::namedWindow(caption, CV_WINDOW_AUTOSIZE);
    cv::imshow(caption, image);
}

inline
int wait_for_key(int msecs)
{
    return cv::waitKey(msecs);
}

// Convert and normalize to double a given cv::Mat image. If the supposed size 
// of pixel (T) differs from real size in the given image, return empty RawImage.
template <typename ValType> template <typename T>
RawImage<ValType> RawImage<ValType>::from_cvmat(const cv::Mat& image)
{
    PixelMatrix retvalue;

    if (is_bpps_equal_to<T>(image))
    {   
        retvalue.resize(image.rows);
        T factor = std::numeric_limits<T>::max();

        for(int i = 0; i < image.rows; ++i)
        {
            retvalue[i].resize(image.cols);

            for(int j = 0; j < image.cols; ++j)
                retvalue[i][j] = static_cast<ValType>(image.at<T>(i, j)) / 
                                     static_cast<ValType>(factor);
        }
    }

    return RawImage(retvalue);
}

// Convert current state to cv::Mat image. CV_8UC1 flag is used to create an image
// with 1 byte per pixel intensities.
template <typename ValType> 
cv::Mat RawImage<ValType>::to_cvmat() const
{
    cv::Mat retvalue(static_cast<int>(image_.size()), 
                     static_cast<int>(image_[0].size()), 
                     CV_8UC1);

    boost::uint8_t factor = std::numeric_limits<boost::uint8_t>::max();

    for(int i = 0; i < retvalue.rows; ++i)
    {
        for(int j = 0; j < retvalue.cols; ++j)
            retvalue.at<boost::uint8_t>(i, j) = static_cast<boost::uint8_t>
                                                    (image_[i][j] * factor);
    }

    return retvalue;
}

#endif


// RawImage methods definition.
//template <typename ValType>
//RawImage<ValType>::RawImage(const PixelMatrix& image):
//    image_(image)
//{ }

//template <typename ValType> inline
//typename RawImage<ValType>::PixelMatrix RawImage<ValType>::raw() const
//{
//    return image_;
//}

template <typename ValType> inline
ValType* RawImage<ValType>::data() const
{
    return image_.get();
}

template <typename ValType> inline
std::size_t RawImage<ValType>::size() const
{
    return (width_ * height_ * sizeof(ValType));
}

//template <typename ValType> inline
//typename const RawImage<ValType>::Pixels& RawImage<ValType>::operator[](std::size_t row) const
//{
//    return image_[row];
//}

//template <typename ValType> inline
//typename RawImage<ValType>::Pixels& RawImage<ValType>::operator[](std::size_t row)
//{
//    return image_[row];
//}

template <typename ValType> inline
ValType& RawImage<ValType>::at(std::size_t col, std::size_t row)
{
    // TODO: Perform checks.
    return image_[col + width_ * row];
}

// Return a set of brightness values of the pixel itself and surrounding neighbours.
template <typename ValType>
typename RawImage<ValType>::Pixels RawImage<ValType>::get_neighbour_values(
    std::size_t col, std::size_t row) const
{
    // TODO: Check index range.

    Pixels retvalue;

    Indices indices = get_neighbours(col, row);
    indices.push_back(std::make_pair(col, row));

    // TODO: Use accessor without check.
    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
        retvalue.push_back(image_.at(it->first, it->second);

    return retvalue;
}

// Returns indices of all first-order neighbours of given pixel.
template <typename ValType>
Indices RawImage<ValType>::get_neighbours(
    std::size_t col, std::size_t row) const
{
    Indices retvalue;

    if (col != 0)
        retvalue.push_back(std::make_pair(col - 1, row));

    if (row != 0)
        retvalue.push_back(std::make_pair(col, row - 1));

    if (col != width_ - 1)
        retvalue.push_back(std::make_pair(col + 1, row));

    if (row != height_ - 1)
        retvalue.push_back(std::make_pair(col, row + 1));

    return retvalue;
}

template <typename ValType>
double RawImage<ValType>::av_dist(std::size_t col, std::size_t row) const
{
    double retvalue = 0.0;

    Indices indices = get_neighbours(col, row);
    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        retvalue += abs(at(col, row) - at(it->first, it->second));
    }

    return
        (retvalue / indices.size());
}

template <typename ValType>
double RawImage<ValType>::std_devia(std::size_t col, std::size_t row) const
{
    double retvalue = 0.0;

    Indices indices = get_neighbours(col, row);

    // First we should get series of differences.
    Pixels diffs;
    diffs.reserve(indices.size());
    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        diffs.push_back(image_[col][row] - image_[it->first][it->second]);
    }

    // Then obtain a mean.
    ValType mean = 0.0;
    for (Pixels::const_iterator it = diffs.begin(); it != diffs.end(); ++it)
    {
        mean += (*it);
    }
    mean = mean / diffs.size();

    // Finally, compute a variance.
    for (Pixels::const_iterator it = diffs.begin(); it != diffs.end(); ++it)
    {
        retvalue += pow((*it) - mean, 2.0);
    }    
    retvalue = sqrt(retvalue) / (diffs.size() - 1);

    return retvalue;
}

} // namespace common
#endif // RAW_IMAGE_HPP_A9C93511_7D52_457E_9B7A_5CFA9590A8C9_
