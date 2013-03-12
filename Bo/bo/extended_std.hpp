
/******************************************************************************

  Extension of the STL I/O streaming.

  Copyright (c) 2012
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

#ifndef EXTENDED_STD_HPP_3122E01D_B758_4638_911C_3F2EF44FC364_
#define EXTENDED_STD_HPP_3122E01D_B758_4638_911C_3F2EF44FC364_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/format.hpp>

namespace bo {

// Streams std::vector<T> contents into a std::string using T::operator<<.
template <typename T>
std::string str(const std::vector<T>& obj)
{
    std::stringstream oss;
    oss << obj;

    return oss.str();
}

} // namespace bo


namespace std {

// Implementation of stream operator<< for std::vector<T>. T must support operator<<.
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& obj)
{
    std::size_t size = obj.size();
    os << boost::format("std::vector of size %1%, object %2$#x: ")
            % size % &obj << std::endl << "    (";

    for (std::size_t i = 0; i + 1 < size; ++i)
        os << boost::format("%1%, %|4t|") % obj[i];

    // Print last element (if any) separately in order to avoid last comma and spaces.
    (obj.empty() ? (os << " )") : (os << boost::format("%1%)") % obj.back())) << std::endl
        << boost::format("end of object %1$#x.") % &obj << std::endl;

    return os;
}

} // namespace std

#endif // EXTENDED_STD_HPP_3122E01D_B758_4638_911C_3F2EF44FC364_
