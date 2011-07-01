
#include <cmath>
#include <gtest/gtest.h>

// Temporary disable warning on conversion (from __int64 to long). This will be fixed
// in the next versions of boost.
#ifdef _MSC_VER
#   pragma warning(push)
#   pragma warning(disable:4244)
#endif // _MSC_VER
#include <boost/thread.hpp>
#ifdef _MSC_VER
#   pragma warning(pop)
#endif // _MSC_VER

#define USE_BOOST_TIMER
#include "common/performance.hpp"

using namespace common;


TEST(TimerTest, BoostTimer)
{
    unsigned sleep_time = 1;

    // This will create boost::timer, since USE_BOOST_TIMER was defined.
    Timer boost_timer;

    // In its current implementation causes C4244 warning. It will be fixed in the
    // next versions of boost.
    boost::this_thread::sleep(boost::posix_time::seconds(sleep_time));


    double elapsed = boost_timer.elapsed();

    // We expect that timer will give an error not more than 20% of the sleep time.
    EXPECT_GE(0.2 * double(sleep_time), abs(elapsed - double(sleep_time)));
}

// This test is only possible under MSVC, since the MSVCTimer class is available
// only for MSVC platform.
#ifdef _MSC_VER

TEST(TimerTest, MSVCTimer)
{
    unsigned sleep_time = 1;

    // This will create MSVCTimer.
    detail::MSVCTimer msvc_timer;

    // In its current implementation causes C4244 warning. It will be fixed in the
    // next versions of boost.
    boost::this_thread::sleep(boost::posix_time::seconds(sleep_time));

    double elapsed = msvc_timer.elapsed();

    // We expect that timer will give an error not more than 20% of the sleep time.
    EXPECT_GE(0.2 * double(sleep_time), abs(elapsed - double(sleep_time)));
}

#endif // _MSC_VER
