
/******************************************************************************

    performance.hpp, v 1.0.1 2011.03.14

    Functions for performance evaluation and manipulation based on WindowsAPI.
    Time is calculated through processor ticks and processor frequency. 

    Copyright (c) 2010, 2011, Alexander Rukletsov <rukletsov@gmail.com>
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

#include <windows.h>

#include <boost/cstdint.hpp>


/* Usage example:
 *
 *    boost::int64_t start = common::get_proc_ticks();
 *
 *    < ... do job here ... > 
 *
 *    boost::int64_t end = common::get_proc_ticks();
 *
 *    std::cout << "job took: "  
 *        << static_cast<float>(end - start) / common::get_proc_freq() * 1000 
 *        << " msecs" << std::endl;
 */


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
