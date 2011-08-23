
/******************************************************************************

    raw_image_2d_io.hpp, v 1.0.0 2011.08.23

    I/O for RawImage2D class. OpenCV library can be used for working with
    image files. This file provides necessary convertion and utility functions
    to work with OpenCV types. Be advised that different versions of OpenCV
    can have different interfaces. Therefore OpenCV 2.0, 2.1 or 2.2 are
    recommended.

    No OpenCV dependency is added to this library. It is possible because all
    functions are templates. However that implies linking against OpenCV for
    users of this file. In other words, if this file is used in a project,
    OpenCV library is needed to be added to the project's dependencies.

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

#ifndef RAW_IMAGE_2D_IO_HPP_203C45E7_9EA0_4B4C_AE8C_A5600E517560_
#define RAW_IMAGE_2D_IO_HPP_203C45E7_9EA0_4B4C_AE8C_A5600E517560_

#include <limits>

// Enable OpenCV usage only when explicitly requested.
// Depending on the version of OpenCV you link to, one of the following
// headers should be used. If your OpenCV version is 2.2, define OPENCV_2_2
// symbol before including this file. If your OpenCV version is either 2.0
// or 2.1, include this file without defining OPENCV_2_2 symbol.
#ifdef LIBCOMMON_USE_OPENCV
#   ifdef OPENCV_2_2
#       include <opencv2/opencv.hpp>
#   else
#       include <opencv/cv.h>
#   endif
#endif

#include "raw_image_2d.hpp"

namespace common {
namespace io {

// Convertion functions from and to OpenCV format are available on demand.
#ifdef LIBCOMMON_USE_OPENCV

// Checks if the size of pixel type is the same as in cv::Mat instance.
template <typename T> inline
bool is_bpps_equal_to(const cv::Mat& model)
{
    return
        (sizeof(T) == model.elemSize());
}

inline
void cv_show_image(const std::string& caption, const cv::Mat& image)
{
    cv::namedWindow(caption, CV_WINDOW_AUTOSIZE);
    cv::imshow(caption, image);
}

inline
int cv_wait_for_key(int msecs)
{
    return cv::waitKey(msecs);
}

// Converts and normalizes to double a given cv::Mat image. If the supposed size
// of pixel differs from real size in the given image, returns empty RawImage2D.
template <typename ValType>
common::RawImage2D<ValType> raw_image_2d_from_cvmat<ValType>(const cv::Mat& image)
{
    if (is_bpps_equal_to<ValType>(image))
    {
        common::RawImage2D<ValType> retvalue(image.cols, image.rows);

        ValType factor = std::numeric_limits<ValType>::max();

        for(int i = 0; i < image.rows; ++i)
        {
            for(int j = 0; j < image.cols; ++j)
                retvalue(i, j) = static_cast<ValType>(image.at<ValType>(i, j)) /
                                     static_cast<ValType>(factor);
        }

        return retvalue;
    }
    else
        return common::RawImage2D<ValType>;
}

// Converts current state to cv::Mat image. CV_8UC1 flag is used to create an image
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

} // namespace io
} // namespace common

#endif // RAW_IMAGE_2D_IO_HPP_203C45E7_9EA0_4B4C_AE8C_A5600E517560_
