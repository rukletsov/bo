
#ifndef DEBUG_ALLOC_HPP_
#define DEBUG_ALLOC_HPP_

// Change new operator to one with more info about the leaked block, let malloc
// remember more info about the leaked block.
#if defined(_MSC_VER) && defined(_DEBUG)
#   define _CRTDBG_MAP_ALLOC
#   include <cstdlib>
#   include <crtdbg.h>
#   define new new(_CLIENT_BLOCK, __FILE__, __LINE__)
#endif // defined(_MSC_VER) && defined(_DEBUG)

#endif // DEBUG_ALLOC_HPP_
