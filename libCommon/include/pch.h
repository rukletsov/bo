// pch.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef PCH_H_
#define PCH_H_

// MSVC-specific stuff and precompiled headers. On platform other than MSVC
// this file will be empty and therefore precompiled headers won't be supported.
#ifdef _MSC_VER

#ifndef _WIN32_WINNT            // Allow use of features specific to Windows XP or later.
#   define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other
#endif                          // versions of Windows.
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

// Specify headers which will be stored in a precompiled header. USE_PCH variable can
// be defined by the user during configuration process (with CMake), or explicitly
// before this line. Note that msvc's pch cannot be used in distributed builds.
#ifdef USE_PCH
#   include <cstddef>
#   include <cmath>
#   include <stdexcept>
#   include <string>
#   include <vector>
#   include <map>
#   include <boost/cstdint.hpp>
#   include <boost/assert.hpp>
#   include <boost/array.hpp>
#   include <boost/format.hpp>
#endif // USE_PCH

#endif // _MSC_VER

#endif // PCH_H_
