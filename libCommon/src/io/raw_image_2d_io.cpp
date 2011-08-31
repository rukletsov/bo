
#include "pch.h"
#include <fstream>

#include "raw_image_2d_io.hpp"


namespace common {
namespace io {

void save_raw_image_float_to_8bpps(common::RawImage2D<float> image,
                                   const std::string& filename)
{
    if (image.is_null())
        return;

    std::ofstream fs(filename.c_str(), std::fstream::out | std::fstream::binary);
    for (std::size_t row = 0; row < image.height(); ++row)
        for (std::size_t col = 0; col < image.width(); ++col)
            fs << boost::uint8_t(image(col, row) * 256);
}

} // namespace io
} // namespace common
