
#include "pch.h"

#include <boost/assert.hpp>
#include <boost/regex.hpp>

#include "ini_reader.hpp"


namespace common {
namespace io {

void IniReader::read_file(const String& file_name)
{
    Line line;
	String cur_section = settings_.default_section;
	String cur_subsection = settings_.default_subsection;

    // add default section and default subsection, they should always present
    add_subsection(settings_.default_section, settings_.default_subsection);

	std::ifstream fin(file_name.c_str());

    while (fin.good()) 
	{
		std::getline(fin, line);
		boost::trim(line);

		if (is_section(line)) 
		{   
			boost::trim_if(line, boost::is_any_of(settings_.section_left_symbols + 
                settings_.section_right_symbols));

			// append 'default_subsection_' to the line
			// using any of the section delimiters
			if (!is_subsection(line))
				line += settings_.section_delimiter_symbols[0] + settings_.default_subsection;

			String section = extract_section(line);
			if (section != settings_.skipped_section)
			{	
				cur_section = section;
				cur_subsection = extract_subsection(line);
				add_subsection(cur_section, cur_subsection);        // adds both section and sub section if any of them does not exist
			} 
        } 
		else if (is_keyvalue(line))
		{   // process property
			String key = extract_key(line);
			String value = extract_value(line);

			// empty keys are not allowed
            if (key != ini::PROHIBITED_KEY)
				add_keyvalue(cur_section, cur_subsection, key, value);
        }
		else if ( line.empty() || is_comments(line) )
		{	// skip the line
		}
		else
		{	// actually, an error; example " key" or "\t some text \t"
		}
    }
}

String IniReader::get_value(const String& section_name, 
							const String& subsection_name, 
							const String& key_name,
							const String& def_value) const
{
	String retvalue = def_value;

	// Check if we have specified section, subsection and key.
	if (has_pair(section_name, subsection_name, key_name))
	{
		retvalue = sections_.find(section_name)->
				   second.find(subsection_name)->
				   second.find(key_name)->second;
	}

    return retvalue;
}

String IniReader::get_value(const String& section_name, 
							const String& key_name, 
							const String& def_value) const
{
	return get_value(section_name, settings_.default_subsection, key_name, def_value);
}

String IniReader::get_value(const String& key_name, 
							const String& def_value) const
{
	return get_value(settings_.default_section, settings_.default_subsection, key_name, def_value);
}

Strings IniReader::get_section_names() const
{
	Strings retvalue;

	for (SectionsIterator sect_it = sections_.begin(); sect_it != sections_.end(); ++sect_it)
		retvalue.push_back(sect_it->first);

	return retvalue;
}

Strings IniReader::get_section_names_by_pattern(const String& section_name_pattern) const
{
    Strings retvalue;

    boost::regex sections_pattern(section_name_pattern);   	
    for (SectionsIterator sect_it = sections_.begin(); sect_it != sections_.end(); ++sect_it)
    {	
        // If current section name matches the pattern, then add it to result
        if (boost::regex_match(sect_it->first, sections_pattern) == true)
            retvalue.push_back(sect_it->first);
    }

	return retvalue;
}

Strings IniReader::get_subsection_names(const String& section_name) const
{
	Strings retvalue;

	SectionsIterator sect_it = sections_.find(section_name);
	if (sect_it != sections_.end())
	{	
		for (SubSectionsIterator subsect_it = sect_it->second.begin(); 
			subsect_it != sect_it->second.end(); ++subsect_it)
		{
			 retvalue.push_back(subsect_it->first);
		}
	}

	return retvalue;
}

Strings IniReader::get_subsection_names_by_pattern(const String& section_name, 
                                                   const String& subsection_name_pattern) const
{
    Strings retvalue;

    boost::regex subsections_pattern(subsection_name_pattern);   	
	SectionsIterator sect_it = sections_.find(section_name);
	if (sect_it != sections_.end())
	{	
		for (SubSectionsIterator subsect_it = sect_it->second.begin(); 
			subsect_it != sect_it->second.end(); ++subsect_it)
		{
			 if (boost::regex_match(subsect_it->first, subsections_pattern) == true)
                 retvalue.push_back(subsect_it->first);
		}
	}

	return retvalue;
}

Strings IniReader::get_key_names(const String &section_name, const String &subsection_name) const
{
	Strings retvalue;

	SectionsIterator sect_it = sections_.find(section_name);
	if (sect_it != sections_.end())
	{	
		SubSectionsIterator subsect_it = sect_it->second.find(subsection_name);
		if (subsect_it != sect_it->second.end())
		{
			for (KeyValuePairsIterator pair_it = subsect_it->second.begin(); 
				pair_it != subsect_it->second.end(); ++pair_it)
			{
				 retvalue.push_back(pair_it->first);
			}
		}
	}

	return retvalue;
}

Strings IniReader::get_key_names_by_pattern(const String& section_name, 
                                            const String& subsection_name,
                                            const String& key_name_pattern) const
{
    Strings retvalue;

    boost::regex keys_pattern(key_name_pattern);   	
	SectionsIterator sect_it = sections_.find(section_name);
	if (sect_it != sections_.end())
	{	
		SubSectionsIterator subsect_it = sect_it->second.find(subsection_name);
		if (subsect_it != sect_it->second.end())
		{
			for (KeyValuePairsIterator pair_it = subsect_it->second.begin(); 
				pair_it != subsect_it->second.end(); ++pair_it)
			{
				if (boost::regex_match(pair_it->first, keys_pattern) == true) 
                    retvalue.push_back(pair_it->first);
			}
		}
	}

	return retvalue;
}

Strings IniReader::get_key_names() const
{
    return 
        get_key_names(settings_.default_section, settings_.default_subsection);
}

Strings IniReader::get_key_names_by_pattern(const String& key_name_pattern) const
{
    return
        get_key_names_by_pattern(settings_.default_section, settings_.default_subsection,
            key_name_pattern);
}

String IniReader::get_default_section_name()
{
    return settings_.default_section;
}

 
String IniReader::get_default_subsection_name()
{
    return settings_.default_subsection;
}

void IniReader::add_subsection(const String& section_name, const String& subsection_name)
{
	// Thanks to std::map.operator[] functionality we can add non-existent 
	// section and/or subsection like this.
	sections_[section_name][subsection_name];
}

void IniReader::add_keyvalue(const String& section_name, 
							 const String& subsection_name, 
							 const String& key_name,
							 const String& value)
{
	// 'section_name' section and 'subsection_name' subsection should already exist in 'sections_'
    BOOST_ASSERT(sections_.find(section_name) != sections_.end() && 
        "IniReader: can not add key/value pair to a non-existent section.");
    BOOST_ASSERT(sections_[section_name].find(subsection_name) != sections_[section_name].end() && 
        "IniReader: can not add key/value pair to a non-existent subsection.");

	sections_[section_name][subsection_name][key_name] = value;
}

Strings IniReader::split_subsection(const Line& trimmed_line) const
{
	Strings retvalue;
	Strings split_result;

	boost::split(split_result, 
				 trimmed_line, 
				 boost::is_any_of(settings_.section_delimiter_symbols), 
				 boost::token_compress_off);

    BOOST_ASSERT(split_result.size() > 1 && "IniReader: subsection name or format is invalid.");

	switch(settings_.section_error_type)
	{
    case IniReaderSettings::TAKE_1ST_AND_SECOND:
		retvalue.push_back(split_result[0]);
		retvalue.push_back(split_result[1]);
		break;

	case IniReaderSettings::TAKE_1ST_AND_LAST:
		retvalue.push_back(split_result.front());
		retvalue.push_back(split_result.back());
		break;

	case IniReaderSettings::IGNORE_SECTION:
		retvalue.push_back(settings_.ignored_section);
		retvalue.push_back(trimmed_line);
		break;

	case IniReaderSettings::IGNORE_LINE:
		retvalue.push_back(settings_.skipped_section);
		retvalue.push_back(settings_.skipped_section);
		break;

	default:
		// this code should be inaccessible
		retvalue.push_back(settings_.erroneous_section);
		retvalue.push_back(trimmed_line);
	}

	return retvalue;
}

Strings IniReader::split_keyvalue(const Line& trimmed_line) const
{
	Strings retvalue;

	switch(settings_.keyvalue_error_type)
	{
	case IniReaderSettings::ALLOW_DELIMITERS_IN_VALUE:
	{
		Line::const_iterator pos = std::find_if(trimmed_line.begin(), 
												trimmed_line.end(), 
												boost::is_any_of(settings_.delimiter_symbols));

		// extract key, do not include delimiter
		String key(trimmed_line.begin(), pos);
		retvalue.push_back(key);

		// extract value, do not include delimiter
		String value(pos + 1, trimmed_line.end());
		retvalue.push_back(value);
		break;
	}
	case IniReaderSettings::SKIP_ERRONEOUS_KEY_VALUE:
	{	
		Strings split_result;
		boost::split(split_result,
			 trimmed_line,
			 boost::is_any_of(settings_.delimiter_symbols),
			 boost::token_compress_off);

        BOOST_ASSERT(split_result.size() > 1 && "IniReader: key/value pair format is invalid.");

		if (split_result.size() > 2)
		{
			retvalue.push_back(ini::PROHIBITED_KEY);
			retvalue.push_back(trimmed_line);
		}
		else
		{	// split_result.size() == 2
			retvalue = split_result;
		}
		break;
	}
	}

	return retvalue;
}

String IniReader::extract_section(const Line& trimmed_line) const
{
	return boost::trim_copy(split_subsection(trimmed_line)[0]);
}

String IniReader::extract_subsection(const Line& trimmed_line) const
{
	return boost::trim_copy(split_subsection(trimmed_line)[1]);
}

String IniReader::extract_key(const Line& trimmed_line) const
{
	return boost::trim_copy(split_keyvalue(trimmed_line)[0]);
}

String IniReader::extract_value(const Line& trimmed_line) const
{
	return boost::trim_copy(split_keyvalue(trimmed_line)[1]);
}

bool IniReader::is_comments(const Line& trimmed_line) const
{
	return
		 (trimmed_line.begin() == std::find_if(trimmed_line.begin(), 
											   trimmed_line.end(), 
											   boost::is_any_of(settings_.comment_symbols)));
}

bool IniReader::is_section(const Line& trimmed_line) const
{
	// Check if the trimmed_line starts with "[".
	return
        boost::starts_with(trimmed_line, settings_.section_left_symbols) && 
        boost::ends_with(trimmed_line, settings_.section_right_symbols);
}

bool IniReader::is_subsection(const Line& trimmed_line) const
{
	// Check if the trimmed_line has one of the section delimiters.
	return 
		(trimmed_line.end() != std::find_if(trimmed_line.begin(), 
											trimmed_line.end(), 
											boost::is_any_of(settings_.section_delimiter_symbols)));
}

bool IniReader::is_keyvalue(const Line& trimmed_line) const
{
	// Check if the trimmed_line has one of the key/value delimiters.
	return
		(trimmed_line.end() != std::find_if(trimmed_line.begin(), 
											trimmed_line.end(), 
											boost::is_any_of(settings_.delimiter_symbols)));		
}

bool IniReader::has_pair(const String& section_name, 
						 const String& subsection_name, 
						 const String& key_name) const
{
	bool retvalue = false;

	SectionsIterator sect_it = sections_.find(section_name);
	if (sect_it != sections_.end())
	{
		SubSectionsIterator subsect_it = sect_it->second.find(subsection_name);
		if (subsect_it != sect_it->second.end())
		{
			KeyValuePairsIterator keyvalue_it = subsect_it->second.find(key_name);
			if (keyvalue_it != subsect_it->second.end())
			{	
				retvalue = true;
			}
		}
	}

	return retvalue;
}

#ifdef _DEBUG
void IniReader::print_sections_() const
{
	std::cout << std::endl << "IniReader dump:" << std::endl;
	for (SectionsIterator sect = sections_.begin(); sect != sections_.end(); ++sect)
	{
		std::cout << sect->first << std::endl;

		for (SubSectionsIterator subsect = sect->second.begin(); subsect != sect->second.end(); ++subsect)
		{
			std::cout << "--" << subsect->first << std::endl;

			for (KeyValuePairsIterator pair = subsect->second.begin(); pair != subsect->second.end(); ++pair)
			{
				std::cout << "----" << pair->first << " = " << pair->second << std::endl;
			}
		}
	}
}
#endif

} // namespace io
} // namespace common
