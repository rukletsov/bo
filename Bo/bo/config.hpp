
/******************************************************************************

  config.hpp, v 1.0.0 2012.03.17

  Library configuration and declarations for porting Bo to various platforms.

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

#ifndef CONFIG_HPP_FE9AAB4E_25DD_4157_A536_7D3DBDD812C3_
#define CONFIG_HPP_FE9AAB4E_25DD_4157_A536_7D3DBDD812C3_

#include <boost/config.hpp>

// Shared library support (dll on Windows).

#ifdef BO_SHARED_LIBRARY
// Export if this is our own source, otherwise import.
#   ifdef BO_SOURCE
#       define BO_DECL BOOST_SYMBOL_EXPORT
#   else
#       define BO_DECL BOOST_SYMBOL_IMPORT
#   endif  // BO_SOURCE
#endif  // BO_SHARED_LIBRARY

#ifndef BO_DECL
#   define BO_DECL
#endif

#endif // CONFIG_HPP_FE9AAB4E_25DD_4157_A536_7D3DBDD812C3_
