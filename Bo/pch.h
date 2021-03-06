
#ifndef BO_PCH
#define BO_PCH

// MSVC-specific stuff.
#ifdef _MSC_VER

// Allow use of features specific to Windows XP or later. Change this to the
// appropriate value to target other versions of Windows. Note that this option can
// be set during configuration process (with CMake). For more information see
//     http://msdn.microsoft.com/en-us/library/aa383745%28v=vs.85%29.aspx
#  ifndef _WIN32_WINNT
#    define _WIN32_WINNT 0x0501
#  endif // _WIN32_WINNT

// Exclude rarely-used stuff from Windows headers. For more information see
//     http://msdn.microsoft.com/en-us/library/aa383745%28v=vs.85%29.aspx
#  define WIN32_LEAN_AND_MEAN

#endif // _MSC_VER

// Specify headers which will be stored in a precompiled header. BO_USE_PCH macro can
// be defined by the user during configuration process (with CMake), or explicitly
// before this line. Note that msvc's pch cannot be used in distributed builds.
// Note that C++ includes should be separated from C includes in order to let the
// precompiled header be used also for C source files.
#ifdef BO_USE_PCH
#  include <limits.h>
#  include <stdlib.h>
#  include <stddef.h>
#  include <math.h>

#  ifdef __cplusplus
#    include <cstddef>
#    include <cmath>
#    include <stdexcept>
#    include <utility>
#    include <string>
#    include <vector>
#    include <map>
#    include <set>
#    include <list>
#    include <algorithm>
#    include <limits>
#    include <functional>
#    include <boost/cstdint.hpp>
#    include <boost/assert.hpp>
#    include <boost/array.hpp>
#    include <boost/format.hpp>
#    include <boost/function.hpp>
#    include <boost/scoped_ptr.hpp>
#    include <boost/operators.hpp>
#    include <boost/regex.hpp>
#    include <boost/scope_exit.hpp>
#  endif // __cplusplus

#endif // BO_USE_PCH

#endif // BO_PCH
