
/******************************************************************************

  logging.hpp, v 1.0.3 2012.09.11

  Routines for logging messages and errors.

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
#include <string>
#include <boost/date_time.hpp>

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

} // namespace bo

#endif // LOGGING_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_
