#ifndef ROUTINES_HPP_5E09DEDF_340B_4A8f_A106_FF79F8E4D888_
#define ROUTINES_HPP_5E09DEDF_340B_4A8f_A106_FF79F8E4D888_

#include <iostream>
#include <string>
#include <windows.h>
#include <math.h>

#include <boost/cstdint.hpp>
#include <opencv/cv.h>
#include <opencv/highgui.h>


namespace common {

inline
void errprint(const std::string& app_name, const std::string& msg)
{
	std::cout << "Error: " << msg << std::endl << std::endl
			  << "Use \"" << app_name << " -h\" for help" << std::endl;
}

inline
boost::int64_t get_proc_ticks()
{
    LARGE_INTEGER retvalue;
    QueryPerformanceCounter(&retvalue);

    return 
        static_cast<boost::int64_t>(retvalue.QuadPart);
}

inline 
boost::int64_t get_proc_freq()
{
    LARGE_INTEGER retvalue;
    QueryPerformanceFrequency(&retvalue);

    return 
        static_cast<boost::int64_t>(retvalue.QuadPart);
}

inline
void increase_priority()
{
    SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
}

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

inline
int round(double x)
{
    return
        static_cast<int>(floor(x + 0.5));
}

} // namespace common

#endif // ROUTINES_HPP_5E09DEDF_340B_4A8f_A106_FF79F8E4D888_
