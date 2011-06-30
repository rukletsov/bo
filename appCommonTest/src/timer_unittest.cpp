
#include <cmath>
#include <boost/thread.hpp>
#include <gtest/gtest.h>

#define USE_BOOST_TIMER
#include "common/performance.hpp"

using namespace common;


TEST(TimerTest, BoostTimer)
{
    unsigned sleep_time = 1000;

    // This will create boost::timer, since USE_BOOST_TIMER was defined.
    Timer boost_timer;
    boost::this_thread::sleep(boost::posix_time::milliseconds(sleep_time));
    double elapsed = boost_timer.elapsed();

    // We expect that timer will give an error not more than 30% of the sleep time.
    EXPECT_GE(0.3 * double(sleep_time), abs(elapsed - double(sleep_time)));
}

// This test is only possible under MSVC, since the MSVCTimer class is available
// only for MSVC platform.
#ifdef _MSC_VER

TEST(TimerTest, MSVCTimer)
{
    unsigned sleep_time = 1000;

    // This will create MSVCTimer.
    detail::MSVCTimer msvc_timer;
    boost::this_thread::sleep(boost::posix_time::milliseconds(sleep_time));
    double elapsed = msvc_timer.elapsed();

    // We expect that timer will give an error not more than 30% of the sleep time.
    EXPECT_GE(0.3 * double(sleep_time), abs(elapsed - double(sleep_time)));
}

#endif // _MSC_VER
