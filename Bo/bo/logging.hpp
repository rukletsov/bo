
/******************************************************************************

  logging.hpp, v 1.0.5 2012.09.12

  Routines for logging classes, messages and errors.

  Copyright (c) 2010 - 2012
  Alexander Rukletsov <rukletsov@gmail.com>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
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

#ifndef LOGGING_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_
#define LOGGING_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/format.hpp>
// Suppress C4127 warning under MSVC while including boost date_time headers.
#ifdef _MSC_VER
#   pragma warning (push)
#   pragma warning (disable:4127)
#   include <boost/date_time.hpp>
#   pragma warning(pop)
#endif // _MSC_VER

namespace bo {

inline
void errprint(const std::string& app_name, const std::string& msg)
{
    std::cout << "Error: " << msg << std::endl << std::endl
              << "Use \"" << app_name << " -h\" for help" << std::endl;
}

inline
void logprint(const std::string& msg)
{
    std::cout << "[" << boost::posix_time::second_clock::local_time().time_of_day()
              << "] " << msg << std::endl;
    std::cout.flush();
}

// Streams std::vector<T> contents into a std::string using T::operator<<.
template <typename T>
std::string str(const std::vector<T>& obj)
{
    std::stringstream oss;

    std::size_t size = obj.size();
    oss << boost::format("std::vector of size %1%, object %2$#x: ")
            % size % &obj << std::endl << "    (";

    for (std::size_t i = 0; i < size - 1; ++i)
        oss << boost::format("%1%, %|4t|") % obj[i];

    // Print last element separately in order to avoid last comma and spaces.
    oss << boost::format("%1%)") % obj[size - 1] << std::endl
        << boost::format("end of object %1$#x.") % &obj << std::endl;

    return oss.str();
}

} // namespace bo

#endif // LOGGING_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_
