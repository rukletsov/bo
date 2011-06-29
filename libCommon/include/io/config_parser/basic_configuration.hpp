
/******************************************************************************

  basic_configuration.hpp, v 1.0.0 2011.05.13

  Basic class for reading and storing data from a configuration file.

  Copyright (c) 2011	
  Alena Bakulina <alena.bakulina@ziti.uni-heidelberg.de>
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


#ifndef BASIC_CONFIGURATION_HPP_5F739AB0_F8FE_4996_ABEB_4E3335AA84C8_
#define BASIC_CONFIGURATION_HPP_5F739AB0_F8FE_4996_ABEB_4E3335AA84C8_

#include "ini_reader.hpp"

namespace common {
namespace io {
namespace config_parser {

class BasicConfiguration
{
public:
    BasicConfiguration() : ini_reader_(DEFAULT_INI_SETTINGS)
    { }

    BasicConfiguration(const IniReaderSettings& settings) : ini_reader_(settings)
    { }

    virtual ~BasicConfiguration() 
    { }

    // Reads a config file using the IniReader class.
	virtual void read_file(const String& file_name);

protected:
    // An object containing data from a config file as a plane text.
    IniReader ini_reader_;

}; // class BasicConfiguration

inline
void BasicConfiguration::read_file(const String& file_name)
{
    ini_reader_.read_file(file_name);
}

} // namespace config_parser
} // namespace io
} // namespace common

#endif // BASIC_CONFIGURATION_HPP_5F739AB0_F8FE_4996_ABEB_4E3335AA84C8_
