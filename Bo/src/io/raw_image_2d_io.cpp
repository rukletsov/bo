
#include "pch.h"
#include <fstream>

#include "bo/io/raw_image_2d_io.hpp"

namespace bo {
namespace io {

void save_raw_image_float_to_8bpps(bo::RawImage2D<float> image,
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

bo::RawImage2D<float> load_raw_image_float_8bpps(const std::string& filename,
                                                 std::size_t width,
                                                 std::size_t height)
{
    // Allocate empty image of requested size.
    bo::RawImage2D<float> image(width, height);

    std::ifstream fs(filename.c_str(), std::fstream::in | std::fstream::binary);
    for (std::size_t row = 0; row < height; ++row)
        for (std::size_t col = 0; col < width; ++col)
            image(col, row) = boost::uint8_t(fs.get() / std::numeric_limits<boost::uint8_t>::max());

    return image;
}

} // namespace io
} // namespace bo
