#ifndef RAW_IMAGE_HPP_563F526B_DA89_4252_9B8C_C26F9F66457C_
#define RAW_IMAGE_HPP_563F526B_DA89_4252_9B8C_C26F9F66457C_

#include <vector>
#include <math.h>

#include <boost/cstdint.hpp>
#include <opencv/cv.h>

#include "routines.hpp"


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

    template <typename T> static
    RawImage<ValType> from_cvmat(const cv::Mat& image);

    cv::Mat to_cvmat() const;
    PixelMatrix raw() const;

    size_t size() const;
    Pixels& operator[](size_t row);
    const Pixels& operator[](size_t row) const;
    ValType& pixval_at(size_t row, size_t col);

    Pixels get_neighbour_values(size_t row, size_t col) const;
    Indices get_neighbours(size_t row, size_t col) const;

protected:
    double av_dist(size_t row, size_t col) const;
    double std_devia(size_t row, size_t col) const;

protected:
    PixelMatrix image_;
};


template <typename ValType> 
RawImage<ValType>::RawImage(const PixelMatrix& image): 
    image_(image)
{ }

// Converts and normalizes to double a given cv::Mat image. If the supposed size 
// of pixel (T) differs from real size in the given image, returns empty RawImage.
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

// Converts current state to cv::Mat image. Uses CV_8UC1 flag to create an image
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

template <typename ValType> 
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
ValType& RawImage<ValType>::pixval_at(size_t row, size_t col)
{
    return image_[row][col];
}


// Returns a set of brightness values of the pixel itself and surrounding neighbours.
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
#endif // RAW_IMAGE_HPP_563F526B_DA89_4252_9B8C_C26F9F66457C_
