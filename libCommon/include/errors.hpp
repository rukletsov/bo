#ifndef ERRORS_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_
#define ERRORS_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_

#include <iostream>
#include <string>


namespace common {

inline
void errprint(const std::string& app_name, const std::string& msg)
{
	std::cout << "Error: " << msg << std::endl << std::endl
			  << "Use \"" << app_name << " -h\" for help" << std::endl;
}

} // namespace common

#endif // ERRORS_HPP_39EEE9C8_E33D_4FF4_9D06_6672AAF8B295_
