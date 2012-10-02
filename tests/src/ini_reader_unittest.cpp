
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "bo/io/ini_reader.hpp"

extern std::string DataDirectory;
static const std::string ini_filename = "ini_reader_test.ini";

using namespace bo::io;


// Create a so-called "text fixture" using base class form GTEST.
// IniReaderTest requires an input ini file, if the file doesn't exist or a file path is
// not defined, the first FileExists test fails and other tests are not run.
class IniReaderTest : public testing::Test
{
protected:
    virtual void SetUp()
    {
        // Get data directory name. Input directory can have or not have trail slashes.
        ini_filepath = boost::filesystem3::path(DataDirectory) /= ini_filename;
    }

    ::testing::AssertionResult IsTestFileAvailable(boost::filesystem3::path filepath)
    {
        return
            (boost::filesystem3::exists(boost::filesystem3::path(filepath))
             ? ::testing::AssertionSuccess() << "file \"" << filepath.string() << "\" found"
             : ::testing::AssertionFailure() << "file \"" << filepath.string() << "\" not found");
    }

    boost::filesystem3::path ini_filepath;
};

TEST_F(IniReaderTest, DefaultConstructor)
{
    IniReader reader;
    EXPECT_STREQ("", reader.get_value("", "", "", "").c_str());
    EXPECT_STREQ("default_value", reader.get_value("some_section", "some_subsection",
                                                   "", "default_value").c_str());
}

TEST_F(IniReaderTest, DefaultBehaviour)
{
    // Check whether the file exists.
    ASSERT_TRUE(IsTestFileAvailable(ini_filepath));

    IniReader reader;
    reader.read_file(ini_filepath.string());

    // Default section, default subsection.
    EXPECT_STREQ("value_no_sec_no_subsec", reader.get_value("key_no_sec_no_subsec", "").c_str());

    // Keys and values with spaces.
    EXPECT_STREQ("value_without_spaces", reader.get_value("key_with_trim_spaces_01", "").c_str());
    EXPECT_STREQ("value_with_trim_spaces", reader.get_value("key_with_trim_spaces_02", "").c_str());

    // Keys with inside spaces.
    EXPECT_STREQ("value_01", reader.get_value("key with inside spaces 01", "").c_str());
    EXPECT_STREQ("value with inside spaces", reader.get_value("key with inside spaces 02", "").c_str());

    // Duplicate keys within one subsection.
    std::string key_duplicate_value = reader.get_value("key_duplicate", "");
    EXPECT_TRUE(("value_duplicate_01" == key_duplicate_value) ||
                ("value_duplicate_02" == key_duplicate_value));

    // Not full key/value pair.
    EXPECT_STREQ("", reader.get_value("key_with_no_value", "").c_str());

    // Different value types.
    EXPECT_EQ(100, reader.get_value<int>("key_int", 0));
    EXPECT_FLOAT_EQ(100.001f, reader.get_value<float>("key_float", 0.0f));
    EXPECT_DOUBLE_EQ(10000.00000100002, reader.get_value<double>("key_double", 0.0f));
    EXPECT_EQ(true, reader.get_value<bool>("key_boolean", 0));

    // Delimiters in value.
    EXPECT_STREQ("value_with_delimiter = second_value_part",
        reader.get_value("key_value_with_delimiters", "").c_str());

    // Case sensitive keys.
    EXPECT_STREQ("value_sensitive_case_small", reader.get_value("key_sensitive_case", "").c_str());
    EXPECT_STREQ("value_sensitive_case_big", reader.get_value("KEY_SENSITIVE_CASE", "").c_str());

    // Section, no subsection.
    EXPECT_STREQ("value_S1", reader.get_value("Section01", "key_S1", "").c_str());

    // Section with one subsection.
    EXPECT_STREQ("value_S2s1",
        reader.get_value("Section02", "Subsection01", "key_S2s1", "").c_str());

    // Section with default subsection and two sections.
    EXPECT_STREQ("value_S3s0",
        reader.get_value("Section03", "key_S3s0", "").c_str());
    EXPECT_STREQ("value_S3s1",
        reader.get_value("Section03", "Subsection01", "key_S3s1", "").c_str());
    EXPECT_STREQ("value_S3s2",
        reader.get_value("Section03", "Subsection02", "key_S3s2", "").c_str());

    // Test with an unsupported subsubsection.
    EXPECT_STREQ("",
        reader.get_value("Section02", "Subsection01.Subsubsection",
                         "key_in_sub_sub_section", "").c_str());
    EXPECT_STREQ("",
        reader.get_value("Section02", "Subsection01", "key_in_sub_sub_section", "").c_str());
    EXPECT_STREQ("value_in_sub_sub_section",
        reader.get_value("Section02", "Subsubsection", "key_in_sub_sub_section", "").c_str());

    // Section with no section.
    EXPECT_STREQ("value_in_subsection_no_section",
        reader.get_value("", "SubsectionNoSection", "key_in_subsection_no_section", "").c_str());

    // Section with a point in the name.
    EXPECT_STREQ("",
        reader.get_value("Section.WithPointInName", "key_in_section_with_a_point_in_name",
                         "").c_str());
    EXPECT_STREQ("value_in_a_section_with_a_point_in_name",
        reader.get_value("Section", "WithPointInName", "key_in_section_with_a_point_in_name",
                         "").c_str());

    // Section with an empty name.
    EXPECT_STREQ("value_in_section_with_an_empty_name",
        reader.get_value("", "key_in_section_with_an_empty_name", "").c_str());

    // Section name with spaces.
    EXPECT_STREQ("value_in_section_with_spaces_in_name",
        reader.get_value("Section name with spaces", "key_in_section_with_spaces_in_name", "").c_str());

    // Subsection before section.
    EXPECT_STREQ("value_in_subsection_before_section",
        reader.get_value("Section03", "SubSectionBeforeSection",
                         "key_in_subsection_before_section", "").c_str());
    EXPECT_STREQ("value_in_subsection_before_section_def_subsection",
        reader.get_value("Section03", "key_in_subsection_before_section_def_subsection", "").c_str());
}

TEST_F(IniReaderTest, OtherDelimiters)
{
    // Check whether the file exists.
    ASSERT_TRUE(IsTestFileAvailable(ini_filepath));

    IniReaderSettings settings;
    settings.comment_symbols = "~";
    settings.delimiter_symbols = "=-";
    settings.section_delimiter_symbols = "_,.";
    settings.section_left_symbols = "<";
    settings.section_right_symbols = ">";

    IniReader reader(settings);
    reader.read_file(ini_filepath.string());

    EXPECT_STREQ("default_value",
        reader.get_value("should_be,ignored-it<is>a~comment_with-delimiters", "default_value").c_str());

    EXPECT_STREQ("value01", reader.get_value("Section01", "key01", "").c_str());

    EXPECT_STREQ("value02_with_section_delimiter",
        reader.get_value("Section01", "Subsection01", "key02", "").c_str());

    EXPECT_FLOAT_EQ(20.5f, reader.get_value<float>("Section02", "Subsection01", "key03", 0.0f));

    EXPECT_STREQ("value04_with<section>delimiters",
        reader.get_value("Section02", "Subsection02", "key04", "").c_str());
}

TEST_F(IniReaderTest, SectionErrorHandlingType)
{
    // Check. whether the file exists.
    ASSERT_TRUE(IsTestFileAvailable(ini_filepath));

    IniReaderSettings settings1;
    settings1.section_error_type = IniReaderSettings::TAKE_1ST_AND_SECOND;
    IniReader reader1(settings1);
    reader1.read_file(ini_filepath.string());

    IniReaderSettings settings2;
    // Parameters should be appended to the previous section.
    settings2.section_error_type = IniReaderSettings::IGNORE_LINE;
    IniReader reader2(settings2);
    reader2.read_file(ini_filepath.string());

    IniReaderSettings settings3;
    // Parameters till the next section are ignored.
    settings3.section_error_type = IniReaderSettings::IGNORE_SECTION;
    IniReader reader3(settings3);
    reader3.read_file(ini_filepath.string());

    // TAKE_1ST_AND_LAST is already tested in DefaultBehaviour test.

    // TAKE_1ST_AND_SECOND
    EXPECT_STREQ("value_in_sub_sub_section",
        reader1.get_value("Section02", "Subsection01", "key_in_sub_sub_section", "").c_str());

    // IGNORE_LINE
    EXPECT_STREQ("value_in_sub_sub_section", reader2.get_value("Section03", "Subsection02",
        "key_in_sub_sub_section", "").c_str());

    // IGNORE_SECTION
    EXPECT_STREQ("",
        reader3.get_value("Section02", "Subsection01", "key_in_sub_sub_section", "").c_str());
    EXPECT_STREQ("",
        reader3.get_value("Section02", "Subsubsection", "key_in_sub_sub_section", "").c_str());
    EXPECT_STREQ("value_in_sub_sub_section", reader3.get_value(settings3.ignored_section,
        "Section02.Subsection01.Subsubsection", "key_in_sub_sub_section", "").c_str());
}

TEST_F(IniReaderTest, KeyValueErrorHandlingType)
{
    // Check whether the file exists.
    ASSERT_TRUE(IsTestFileAvailable(ini_filepath));

    // 1. ALLOW_DELIMITERS_IN_VALUE case tested in DefaultBehaviour test.

    // 2. SKIP_ERRONEOUS_KEY_VALUE
    IniReaderSettings settings2;
    settings2.keyvalue_error_type = IniReaderSettings::SKIP_ERRONEOUS_KEY_VALUE;
    IniReader reader2(settings2);
    reader2.read_file(ini_filepath.string());

    EXPECT_STREQ("default_value", reader2.get_value(ini::PROHIBITED_KEY, "default_value").c_str());
}
