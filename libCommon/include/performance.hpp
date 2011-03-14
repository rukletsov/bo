#ifndef PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_
#define PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_

#include <windows.h>

#include <boost/cstdint.hpp>


namespace common {

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
void increase_process_priority()
{
    SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
}

} // namespace common

#endif // PERFORMANCE_HPP_88D42F69_941B_451A_BC38_DAF8A399F977_
