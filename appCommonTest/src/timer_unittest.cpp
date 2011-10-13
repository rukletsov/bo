
#include <cmath>
#include <ctime>
#include <gtest/gtest.h>

// This allows us to have (and to test) both timers under MSVC: boost::timer through
// Timer class and a special high-resolution timer through detail::MSVCTimer.
#define USE_BOOST_TIMER
#include "common/performance.hpp"

#include "debug_alloc.hpp"

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

    // Perform time queries for a giving period of time. Suspending thread using
    // boost::this_thread::sleep() doesn't work on Linux, since std::clock() counts
    // only working time, not sleeping time.
    std::clock_t endwait = std::clock() + sleep_time * CLOCKS_PER_SEC ;
    while (clock() < endwait)
    { }

    double elapsed = timer.elapsed();

    // We expect that timer will give an error not more than 20% of the sleep time.
    EXPECT_GE(0.2 * double(sleep_time), abs(elapsed - double(sleep_time)));
}
