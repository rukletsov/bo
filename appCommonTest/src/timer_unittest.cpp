
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

// This allows us to have (and to test) both timers under MSVC: boost::timer through
// Timer class and a special high-resolution timer through detail::MSVCTimer.
#define USE_BOOST_TIMER
#include "common/performance.hpp"

using namespace common;


// Create a fixture class template for all timer classes.
template <typename TimerType>
class TimerTest: public testing::Test
{ };

// Associate a list of available timers with the test case. Testing of the MSVCTimer
// class is only possible under MSVC, since the MSVCTimer class is available only for
// MSVC platform.
#ifdef _MSC_VER
    typedef testing::Types<Timer, detail::MSVCTimer> TimerTypes;
#else // _MSC_VER
    typedef testing::Types<Timer> TimerTypes;
#endif // _MSC_VER

TYPED_TEST_CASE(TimerTest, TimerTypes);


TYPED_TEST(TimerTest, MeasuringSleep)
{
    unsigned sleep_time = 1;

    // Create an instance of some timer, which type is one of the specified above.
    TypeParam timer;

    // In its current implementation causes C4244 warning. It will be fixed in the
    // next versions of boost.
    boost::this_thread::sleep(boost::posix_time::seconds(sleep_time));

    double elapsed = timer.elapsed();

    // We expect that timer will give an error not more than 20% of the sleep time.
    EXPECT_GE(0.2 * double(sleep_time), abs(elapsed - double(sleep_time)));
}
