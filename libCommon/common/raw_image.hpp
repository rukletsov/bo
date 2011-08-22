
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


typedef std::pair<size_t, size_t> Index;
typedef std::vector<Index> Indices;


template <typename ValType>
class RawImage
{
public:
    typedef std::vector<ValType> Pixels;
    typedef std::vector<Pixels> PixelMatrix;

    RawImage(const PixelMatrix& image);

    // Convertion functions from and to OpenCV format are available on demand.
#ifdef USE_OPENCV
    template <typename T> static
    RawImage<ValType> from_cvmat(const cv::Mat& image);

    cv::Mat to_cvmat() const;
#endif

    PixelMatrix raw() const;

    size_t size() const;
    Pixels& operator[](size_t row);
    const Pixels& operator[](size_t row) const;
    ValType& at(size_t row, size_t col);

    Pixels get_neighbour_values(size_t row, size_t col) const;
    Indices get_neighbours(size_t row, size_t col) const;

protected:
    double av_dist(size_t row, size_t col) const;
    double std_devia(size_t row, size_t col) const;

protected:
    PixelMatrix image_;
};


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
template <typename ValType> 
RawImage<ValType>::RawImage(const PixelMatrix& image): 
    image_(image)
{ }

template <typename ValType> inline
typename RawImage<ValType>::PixelMatrix RawImage<ValType>::raw() const
{
    return image_;
}

template <typename ValType> inline
size_t RawImage<ValType>::size() const
{
    return image_.size();
}

template <typename ValType> inline
typename const RawImage<ValType>::Pixels& RawImage<ValType>::operator[](size_t row) const
{
    return image_[row];
}

template <typename ValType> inline
typename RawImage<ValType>::Pixels& RawImage<ValType>::operator[](size_t row)
{
    return image_[row];
}

template <typename ValType> inline
ValType& RawImage<ValType>::at(size_t row, size_t col)
{
    return image_[row][col];
}


// Return a set of brightness values of the pixel itself and surrounding neighbours.
template <typename ValType>
typename RawImage<ValType>::Pixels RawImage<ValType>::get_neighbour_values(
    size_t row, size_t col) const
{
    Pixels retvalue;

    Indices indices = get_neighbours(row, col);
    indices.push_back(std::make_pair(row, col));

    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
        retvalue.push_back(image_[it->first][it->second]);

    return retvalue;
}

// Returns indices of all first-order neighbours of given pixel.
template <typename ValType>
Indices RawImage<ValType>::get_neighbours(
    size_t row, size_t col) const
{
    Indices retvalue;

    if (row != 0)
        retvalue.push_back(std::make_pair(row - 1, col));

    if (col != 0)
        retvalue.push_back(std::make_pair(row, col - 1));

    if (row != image_.size() - 1)
        retvalue.push_back(std::make_pair(row + 1, col));

    if (col != image_[row].size() - 1)
        retvalue.push_back(std::make_pair(row, col + 1));

    return retvalue;
}

template <typename ValType>
double RawImage<ValType>::av_dist(size_t row, size_t col) const
{
    double retvalue = 0.0;

    Indices indices = get_neighbours(row, col);
    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        retvalue += abs(image_[row][col] - image_[it->first][it->second]);
    }

    return
        (retvalue / indices.size());
}

template <typename ValType>
double RawImage<ValType>::std_devia(size_t row, size_t col) const
{
    double retvalue = 0.0;

    Indices indices = get_neighbours(row, col);

    // First we should get series of differences.
    Pixels diffs;
    diffs.reserve(indices.size());
    for (Indices::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        diffs.push_back(image_[row][col] - image_[it->first][it->second]);
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
