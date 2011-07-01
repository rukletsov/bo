
/******************************************************************************

    performance.hpp, v 1.1.2 2011.06.30

    Timer class for performance evaluation. By default uses boost::timer
    class. However on Windows a special more precise alternative can be used
    (based on performance counter and Windows API).

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

#ifndef PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_
#define PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_

#include <boost/timer.hpp>

#ifdef _MSC_VER
#   include <limits>
#   define NOMINMAX
#   include <windows.h>
#endif // _MSC_VER

/** The interface of the Timer class is identical with the one of boost::timer class.
  * On the platforms other than Windows Timer class is just a typedef for boost::timer.
  *
  * On Windows a more precise technique can be used: performance counter, which is
  * implemented through MSVCTimer class. Its interface is identical with the one of 
  * boost::timer class. On Windows (actually, where _MSC_VER symbol is defined) by
  * default Timer class coincides with MSVCTimer. For more information see
  *     http://msdn.microsoft.com/en-us/library/ms644904%28v=vs.85%29.aspx
  *     http://msdn.microsoft.com/en-us/library/ms644905%28v=vs.85%29.aspx
  *
  * However the installed hardware can lack the support of the high-resolution
  * performance counter. In this case MSVCTimer works incorrectly and function
  * MSVCTimer::is_supported() returns false. A standard boost::timer should be used
  * instead. To supress default behaviour and use boost::timer instead of MSVCTimer 
  * on Windows (where _MSC_VER is defined), define USE_BOOST_TIMER before including
  * this header.
  *
  * Usage example:
  *
  *     Timer timer;
  *
  *     < ... do calculations here ... >
  *
  *     std::cout << "calculations took: "
  *               << timer.elapsed() << " seconds" << std::endl;
  *
  */


namespace common {

#if !defined(_MSC_VER) || defined(USE_BOOST_TIMER)
    // Use boost::timer on non-Windows and by default.
    typedef boost::timer Timer;
#else 
    // Use MSVCTimer otherwise.
    class detail::MSVCTimer;
    typedef detail::MSVCTimer Timer;
#endif // !defined(_MSC_VER) || defined(USE_BOOST_TIMER)


#ifdef _MSC_VER

namespace detail {

class MSVCTimer
{
public:
    MSVCTimer();
    void restart();
    double elapsed() const;
    double elapsed_max() const;
    double elapsed_min() const;

    static bool is_supported();

protected:
    LONGLONG get_proc_ticks_() const;
    LONGLONG get_proc_freq_() const;

protected:
    LONGLONG start_time_;
};


inline
MSVCTimer::MSVCTimer()
{
    restart();
}

inline
void MSVCTimer::restart()
{
    start_time_ = get_proc_ticks_();
}

inline
double MSVCTimer::elapsed() const
{
    return
        (double(get_proc_ticks_() - start_time_)) / get_proc_freq_();
}

inline
double MSVCTimer::elapsed_max() const
{
    return
        (double(std::numeric_limits<LONGLONG>::max() - start_time_)) /
            get_proc_freq_();
}

inline
double MSVCTimer::elapsed_min() const
{
    return
        double(1) / get_proc_freq_();
}

inline
bool MSVCTimer::is_supported()
{
    LARGE_INTEGER temp;
    return
        (0 == QueryPerformanceFrequency(&temp) ? false : true);
}

inline
LONGLONG MSVCTimer::get_proc_ticks_() const
{
    LARGE_INTEGER retvalue;
    QueryPerformanceCounter(&retvalue);

    return retvalue.QuadPart;
}

inline
LONGLONG MSVCTimer::get_proc_freq_() const
{
    LARGE_INTEGER retvalue;
    QueryPerformanceFrequency(&retvalue);

    return retvalue.QuadPart;
}

} // namespace detail

#endif // _MSC_VER

} // namespace common

#endif // PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_
