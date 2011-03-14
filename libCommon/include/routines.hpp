#ifndef ROUTINES_HPP_F5FE5D7E_04CE_4B3C_AC4A_6A8B6D478801_
#define ROUTINES_HPP_F5FE5D7E_04CE_4B3C_AC4A_6A8B6D478801_

#include <math.h>

#include <boost/cstdint.hpp>


namespace common {

inline
int round(double x)
{
    return
        static_cast<int>(floor(x + 0.5));
}

} // namespace common

#endif // ROUTINES_HPP_F5FE5D7E_04CE_4B3C_AC4A_6A8B6D478801_
