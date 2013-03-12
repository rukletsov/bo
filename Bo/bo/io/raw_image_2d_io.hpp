
/******************************************************************************

  I/O for RawImage2D class. OpenCV library can be used for working with
  image files. This file provides necessary convertion and utility functions
  to work with OpenCV types. Be advised that different versions of OpenCV
  can have different interfaces. Therefore OpenCV 2.0, 2.1 or 2.2 are
  recommended.

  No OpenCV dependency is added to this library. It is possible because all
  functions are templates. However that implies linking against OpenCV for
  users of this file. In other words, if this file is used in a project,
  OpenCV library is needed to be added to the project's dependencies.

  Copyright (c) 2010 - 2013
  Alexander Rukletsov <rukletsov@gmail.com>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
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

#include "bo/config.hpp"

#include <limits>
#include <string>
#include <boost/cstdint.hpp>
#include <fstream>

#include "bo/raw_image_2d.hpp"

// Enable OpenCV usage only when explicitly requested.
// Depending on the version of OpenCV you link to, one of the following
// headers should be used. If your OpenCV version is 2.2, define OPENCV_2_2
// symbol before including this file. If your OpenCV version is either 2.0
// or 2.1, include this file without defining OPENCV_2_2 symbol.
#ifdef BO_USE_OPENCV
#   ifdef OPENCV_2_2
#       include <opencv2/opencv.hpp>
#   else
#       include <opencv/cv.h>
#   endif
#endif

namespace bo {
namespace io {

// Saves normalized RawImage2D<RealType> to a raw file (8 bytes per pixel, row by row).
template <typename RealType>
void save_raw_image_to_8bpps(bo::RawImage2D<RealType> image,
                             const std::string& filename)
{
    if (image.is_null())
        return;

    std::ofstream fs(filename.c_str(), std::fstream::out | std::fstream::binary);
    for (std::size_t row = 0; row < image.height(); ++row)
        for (std::size_t col = 0; col < image.width(); ++col)
            fs << boost::uint8_t(image(col, row) *
                                 std::numeric_limits<boost::uint8_t>::max());
}

// Loads normalized RawImage2D<RealType> from a raw file (8 bytes per pixel, row by row).
template <typename RealType>
bo::RawImage2D<RealType> load_raw_image_8bpps(const std::string& filename,
                                              std::size_t width,
                                              std::size_t height)
{
    // Allocate empty image of requested size.
    bo::RawImage2D<RealType> image(width, height);

    std::ifstream fs(filename.c_str(), std::fstream::in | std::fstream::binary);
    for (std::size_t row = 0; row < height; ++row)
        for (std::size_t col = 0; col < width; ++col)
            image(col, row) = RealType(boost::uint8_t(fs.get() /
                                       std::numeric_limits<boost::uint8_t>::max()));

    return image;
}

// Helper functions for OpenCV library are available on demand.
#ifdef BO_USE_OPENCV

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

#endif

} // namespace io
} // namespace bo

#endif // RAW_IMAGE_2D_IO_HPP_203C45E7_9EA0_4B4C_AE8C_A5600E517560_
