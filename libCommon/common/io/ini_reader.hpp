
/******************************************************************************

	ini_reader.hpp, v 1.0.2 2011.05.16

	IniReader implemetation for plain ini-files with common syntax.

    Copyright (c) 2009-2011
    Alexander Rukletsov <rukletsov@gmail.com>
    Alexander Gusak <fami.alex.lom@gmail.com>
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
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODSll
	OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
	HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICTl
	LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
	OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
	SUCH DAMAGE.

*******************************************************************************/


#ifndef INI_READER_HPP_75719B42_4225_4263_8F0C_29EA63E1A2B5_
#define INI_READER_HPP_75719B42_4225_4263_8F0C_29EA63E1A2B5_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/algorithm/string.hpp>

// Used only in print_sections_() function.
#ifdef _DEBUG
    #include <iostream>
#endif

namespace common {
namespace io {

typedef std::string String;
typedef std::string Symbols;
typedef std::string Line;
typedef std::vector<String> Strings;

namespace ini {

// We cannot define not integral constants inside the class. That's why they are here.
static const Symbols COMMENTS = "#;";
static const Symbols DELIMITERS = "=:";
static const Symbols SECTION_DELIMITERS = ".";
static const Symbols SECTION_LEFT = "[";
static const Symbols SECTION_RIGHT = "]";

static const String DEFAULT_SECTION_NAME = "";
static const String DEFAULT_SUBSECTION_NAME = "";
static const String IGNORED_SECTION_NAME = ".ignored";
static const String SKIPPED_SECTION_NAME = ".skipped";
static const String ERRONEOUS_SECTION_NAME = ".erroneous";
static const String PROHIBITED_KEY = "";

} // namespace ini

/* Structure IniReaderSettings provides settings for IniReader class,
 defining it's (reader's) behaviour when parsing file: delimiters 
 symbols, default values, erroneous situations handling types.
 Each delimiter should consist of one symbol (*).
    
 (*) - see comments at IniReader class description.
*/
struct IniReaderSettings
{
public:
    // This enum specifies how to handle errors when subsection name
    // has more than one delimiter, example: "[Section.Sub.Section.1]".
    enum SECTION_ERROR_HANDLING_TYPE
    {
        TAKE_1ST_AND_SECOND,	// --> "Section", "Sub"
        TAKE_1ST_AND_LAST,		// --> "Section", "1"
        IGNORE_SECTION,			// all parameters till the next section are ignored
        IGNORE_LINE				// all parameters will be appended to the previous section
    };

    // This enum specifies how to handle errors when processing key/value line.
    // Example: let key/value delimiter be '=', consider "key = val =ue".
    enum KEY_VALUE_ERROR_HANDLING_TYPE
    {
        ALLOW_DELIMITERS_IN_VALUE,
        SKIP_ERRONEOUS_KEY_VALUE
    };

    // Constructor
    IniReaderSettings(const Symbols& _comment_symbols = ini::COMMENTS, 
			          const Symbols& _delimeter_symbols = ini::DELIMITERS,
			          const Symbols& _section_delimiter_symbols = ini::SECTION_DELIMITERS, 
			          const Symbols& _section_left_symbols = ini::SECTION_LEFT,
			          const Symbols& _section_right_symbols = ini::SECTION_RIGHT,
                      const String& _default_section = ini::DEFAULT_SECTION_NAME,
			          const String& _default_subsection = ini::DEFAULT_SUBSECTION_NAME,
			          const String& _ignored_section = ini::IGNORED_SECTION_NAME,
			          const String& _skipped_section = ini::SKIPPED_SECTION_NAME,
			          const String& _erroneous_section = ini::ERRONEOUS_SECTION_NAME,
                      SECTION_ERROR_HANDLING_TYPE _section_error_type = TAKE_1ST_AND_LAST,
                      KEY_VALUE_ERROR_HANDLING_TYPE _keyvalue_error_type = 
                        ALLOW_DELIMITERS_IN_VALUE):
        comment_symbols(_comment_symbols),
		delimiter_symbols(_delimeter_symbols),
		section_delimiter_symbols(_section_delimiter_symbols),
		section_left_symbols(_section_left_symbols),
		section_right_symbols(_section_right_symbols),
		default_section(_default_section),
        default_subsection(_default_subsection),
		ignored_section(_ignored_section),
		skipped_section(_skipped_section),
		erroneous_section(_erroneous_section),
		section_error_type(_section_error_type),
		keyvalue_error_type(_keyvalue_error_type)
        { }

    // Destructor
    virtual ~IniReaderSettings() { }

public:
    // Symbols used to separate values of different meaning.
	// Every symbols set should contain at least one symbol.
	// Set these fields through constructor to satisfy your requirement.
	Symbols comment_symbols;
	Symbols delimiter_symbols;
	Symbols section_delimiter_symbols;
	Symbols section_left_symbols;
	Symbols section_right_symbols;

    // 'key'/'value' pairs can be places outside any section (in the beginning of
    // a file before any section). Such pairs will be placed in a section with
    // 'default_section' name.
    String default_section;

	// Section can contain 'key'/'value' pairs not nested in a subsection.
	// Such pairs will be placed in subsection with 'default_subsection' name.
	String default_subsection;

	// Ignored sections (if 'section_error_type' is IGNORE_SECTION) will be stored 
	// in the section with this name. If you have changed 
	// 'section_delimiter_symbols' so, that '.' symbol is allowed in 
	// sections' name, change this value in order to avoid possible name collision.
	String ignored_section;

	// Skipped sections (if 'section_error_type' is IGNORE_LINE) are marked with 
	// this name. Be sure no section with this name is possible in an ini-file.
	String skipped_section;

	// All unhandled erroneous sections will be stored here.
	// This section is only for debug. It should never exist in real apps.
	String erroneous_section;

	// Specifies how to handle parsing errors. See comments on 
	// 'SECTION_ERROR_HANDLING_TYPE' type for more information.
	SECTION_ERROR_HANDLING_TYPE section_error_type;

	// Specifies how to handle parsing errors. See comments on 
	// 'KEY_VALUE_ERROR_HANDLING_TYPE' type for more information.
	KEY_VALUE_ERROR_HANDLING_TYPE keyvalue_error_type;
};

// Default settings object.
static const IniReaderSettings DEFAULT_INI_SETTINGS;


/* Class IniReader supports normal ini-file syntax:
 'key'/'value' pairs, comments, blank lines, sections.
 Delimiters and comments' symbols are taken from IniReaderSettings
 class object. Each delimiter should consist of one symbol (*). 
 Sections can contain subsections, but subsections can not. 
 Sections can be omitted, then default section is used.
 Sections and subsections can be declared several times in an ini-file.
 Whitespaces (if any) before 'key', between 'key' and 'delimiter-symbol', 
 between 'delimiter-symbol' and 'value', trailing whitespaces after 
 'value' are ignored. Empty keys are not supported i.e. "  = value" 
 is a bad formed line. Multiple keys inside a subsection are not 
 allowed by now. If an ini-file contains some, one of them will be chosen.

 (*) - You can provide support for multiplied delimiters, such as
	   '==' or ':::' by changing 'boost::token_compress_off' to
	   'boost::token_compress_on' in appropriate split functions.

 TODO:
		1. Add handling for multiple keys (multimap?);
*/

class IniReader
{
public:

    IniReader(const IniReaderSettings& settings = DEFAULT_INI_SETTINGS):
        settings_(settings)
    { }

    virtual ~IniReader() 
    { }

	// Parses ini-file and fills 'sections_' with data.
	void read_file(const String& file_name);

	// Returns 'value' by 'section_name', 'subsection_name', and 'key_name'.
	// If specified section, subsection or key doesn't exist returns 'def_value'.
	// Uses std::stringstream to convert values from String to T.
    template <typename T> T get_value(const String& section_name, 
									  const String& subsection_name, 
									  const String& key_name, 
									  const T& def_value) const;

	// get_value() for 'default_subsection_'.
	template <typename T> T get_value(const String& section_name, 
									  const String& key_name,
									  const T& def_value) const;

    // get_value() for 'default_section_' and 'default_subsection_'.
	template <typename T> T get_value(const String& key_name,
									  const T& def_value) const;

	// Specialization of 4-params get_value() for String values.
	// In our case std::stringstream semantics for String >> String is not suitable.
	String get_value(const String& section_name, 
					 const String& subsection_name, 
					 const String& key_name, 
					 const String& def_value) const;

	// Specialization of 4-params get_value() for String values and
	// when 'subsection_name' is 'default_subsection_'.
	String get_value(const String& section_name, 
					 const String& key_name, 
					 const String& def_value) const;

    // Specialization of 4-params get_value() for String values
	// when 'section_name' is 'default_section_' and
    // 'subsection_name' is 'default_subsection_'.
	String get_value(const String& key_name, 
					 const String& def_value) const;

	// Returns names of all sections, including ignored and erroneous.
	Strings get_section_names() const;

    // Returns names of all sections, including ignored and erroneous,
    // which match an input section name pattern.
    // boost::regex_match is used for matching.
	Strings get_section_names_by_pattern(const String& section_name_pattern) const;

	// Returns names of all subsections incide a given section.
	// If there is no such section returns empty vector.
	Strings get_subsection_names(const String& section_name) const;

    // Returns names of all subsections incide a given section,
	// which match an input subsection name pattern.
	// If there is no such section returns empty vector.
    Strings get_subsection_names_by_pattern(const String& section_name, 
                                            const String& subsection_name_pattern) const;

	// Returns names of all keys incide a given section and subsection.
	// If there is no such section or subsection returns empty vector.
	Strings get_key_names(const String& section_name, 
                          const String& subsection_name) const;

    // Returns names of all keys incide a given section and subsection,
    // which match an input key name pattern.
	// If there is no such section or subsection returns empty vector.
    Strings get_key_names_by_pattern(const String& section_name, 
                                     const String& subsection_name,
                                     const String& key_name_pattern) const;

    // Returns names of all keys incide an default section and default subsection.
	// If there is no such section or subsection returns empty vector.
	Strings get_key_names() const;

    // Returns names of all keys incide an default section and default subsection,
	// which match an input key name pattern.
    // If there is no such section or subsection returns empty vector.
	Strings get_key_names_by_pattern(const String& key_name_pattern) const;

    // Returns a name of a default section.
    String get_default_section_name();

    // Returns a name of a default subsection.
    String get_default_subsection_name();

#ifdef _DEBUG
	void print_sections_() const;
#endif

private:

	typedef std::map<String, String> KeyValuePairs;
	typedef std::map<String, KeyValuePairs> SubSections;
	typedef std::map<String, SubSections> Sections;

	typedef Sections::const_iterator SectionsIterator;
	typedef SubSections::const_iterator SubSectionsIterator;
	typedef KeyValuePairs::const_iterator KeyValuePairsIterator;

	// Adds subsection by section name and subsection name.
	void add_subsection(const String& section_name, const String& subsection_name);

	// Adds key/value pair into a corresponding subsection.
    void add_keyvalue(const String& section_name, 
					  const String& subsection_name, 
					  const String& key_name, 
					  const String& value);

	// Splits subsection name using trimmed line with subsection name from a file.
	// Uses section_error_type_ and other parameters to perform the operation.
	Strings split_subsection(const Line& trimmed_line) const;

	// Splits key/value pair using trimmed line with from a file.
	// Uses 'keyvalue_error_type_' and other parameters to perform the operation.
	Strings split_keyvalue(const Line& trimmed_line) const;

	String extract_section(const Line& trimmed_line) const;

	String extract_subsection(const Line& trimmed_line) const;

	String extract_key(const Line& trimmed_line) const;

	String extract_value(const Line& trimmed_line) const;

	// Checks if the line starts with one of the comment symbols.
	bool is_comments(const Line& trimmed_line) const;

	// Checks if the line contains valid section or subsection name.
	bool is_section(const Line& trimmed_line) const ;

	// Checks if the line contains subsection name (in opposite to section name).
	// Supposes the line is valid section.
	bool is_subsection(const Line& trimmed_line) const;

	// Checks if the line is true key/value pair. Empty key is not checked here.
	bool is_keyvalue(const Line& trimmed_line) const;

	bool has_pair(const String& section_name, const String& subsection_name, 
					  const String& key_name) const;

private:

	// Container for data.
    Sections sections_;

    // Settings defining strings parsing behaviour, delimiters and format.
    IniReaderSettings settings_;
};


template <typename T>
T IniReader::get_value(const String& section_name, 
					   const String& subsection_name, 
					   const String& key_name,
					   const T& def_value) const
{
	T retvalue = def_value;

	// Check if we have specified section, subsection and key.
	if (has_pair(section_name, subsection_name, key_name))
	{	// now we can try to convert 'value' to specified type
		try
		{	
			std::stringstream ss(sections_.find(section_name)->
								 second.find(subsection_name)->
								 second.find(key_name)->second);
			ss >> retvalue;
		}
		catch (...)
		{
			retvalue = def_value;
		}
	}

    return retvalue;
}

template <typename T> inline
T IniReader::get_value(const String& section_name, 
					   const String& key_name,
					   const T& def_value) const
{
	return 
        get_value(section_name, settings_.default_subsection, key_name, def_value);
}

template <typename T> inline
T IniReader::get_value(const String& key_name,
					   const T& def_value) const
{
	return 
        get_value(settings_.default_section, settings_.default_subsection, 
                  key_name, def_value);
}

} // namespace io
} // namespace common

#endif // INI_READER_HPP_75719B42_4225_4263_8F0C_29EA63E1A2B5_
