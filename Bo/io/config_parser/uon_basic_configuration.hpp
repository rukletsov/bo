
/******************************************************************************

  Basic template class for reading and storing data from a config file in 
  the basic UON (UltraOsteon) format (*.conf file). 
  Supports Bo's vector and matrix types. Suppors Atk library types
  if the ATK_SUPPORTED directive is defined.

  Copyright (c) 2011
  Alena Bakulina <alena.bakulina@ziti.uni-heidelberg.de>
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

#ifndef UON_BASIC_CONFIGURATION_HPP_7A47DDA3_838A_4177_AFD0_227E7B76FBAB_
#define UON_BASIC_CONFIGURATION_HPP_7A47DDA3_838A_4177_AFD0_227E7B76FBAB_

#include <vector>
#include <boost/algorithm/string/detail/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "bo/core/vector.hpp"
#include "bo/io/ini_reader.hpp"
#include "bo/io/config_parser/basic_configuration.hpp"
#include "bo/math/blas_extensions.hpp"

// Enable ATK usage only when requested by the library user.
#ifdef ATK_SUPPORTED
#include "utils/vec3.h"
#include "utils/vec4.h"
#include "utils/mat4.h"
#endif

/* Usage example:

    < ... >
    #define ATK_SUPPORTED       // define if you need an Atk data types support
    #include "uon_basic_configuration.hpp"
    < ... >

    < ... some function body ... >

    std::string filename = "..\\_data\\transform0001.conf";
    
    // Config Parser for Bo data types
    bo::io::config_parser::UonBasicConfig uc;
    uc.read_file(filename);
    uc.parse();

    // Read data from the parser
    // See UonBasicConfig definition for other datafields
    std::vector<bo::Vector<float, 3>> points_3d = uc.contour_points_3d;
    math::matrix<float> transform = uc.image_transformation;
    < ... do job with the data here ... >

    // Config Parser for Atk data types
    #ifdef ATK_SUPPORTED
    bo::io::config_parser::UonBasicConfigAtk uc_atk;
    uc_atk.read_file(filename);
    uc_atk.parse();
    // Read data from the parser
    std::vector<vec3f> points_3d_atk = uc_atk.contour_points_3d;
    mat4f transform_atk = uc_atk.image_transformation;
    < ... do job with the data here ... >
    #endif

    < ... >
}
*/

namespace bo {
namespace io {
namespace config_parser {

/* UON basic configuration format:
 UON basic configuration stores vectors and transformation matrices of float
 values. Vectors can be 2D or 3D, but all are stored in fields of a 3D vector type.
 In an ini file a vector is written in a format:
    <key_name> = [f f f] (for a 3D vector, or [f f] for a 2D case).
 Matrices are represented in a file by their's 3D column vectors, each column
 vector has it's own key in a file. Each matrix is defined by four vectors:
 rotational vectors by X, Y and Z and a translation vector.
*/

// Bo's standard matrix type.
typedef math::matrix<float> BoFloatMatrix;


// Namespace uon_config contains definitions of key names for UON configurations.
namespace uon_config {

// Symbols at the beginning and at the end of a value representing a vector (in an
// ini file).
static const Symbols VECTOR_TRIM_SYMBOLS = "[]";

// Names of the keys in a config file.
static const String CONTOUR_POINT_2D_NAME_PATTERN = "contourPoint\\d{4,}2d";
static const String CONTOUR_POINT_3D_NAME_PATTERN = "contourPoint\\d{4,}3d";
static const String IMAGE_TRANSFORM_ROT_X_NAME = "ImageTransform.Rot.x";
static const String IMAGE_TRANSFORM_ROT_Y_NAME = "ImageTransform.Rot.y";
static const String IMAGE_TRANSFORM_ROT_Z_NAME = "ImageTransform.Rot.z";
static const String IMAGE_TRANSFORM_TRANS_NAME = "ImageTransform.Trans";
static const String SUB_IMAGE_IN_PIX_COS_UL_NAME = "SubImageInPixCosUL";
static const String SUB_IMAGE_IN_PIX_COS_UR_NAME = "SubImageInPixCosUR";
static const String IMAGE_SCALING_NAME = "ImageScaling";
static const String WATER_MEAT_COMPENSATION_NAME = "WaterMeatCompensation";

} // namespace uon_config


/* UonBasicConfiguration_<> class is an abstract class which is used as an 
 ancestor of UonBasicConfiguration<> class to provide sharing data fields for
 all specializations of UonBasicConfiguration<> without code repetition.
*/
template <class VectorType, class MatrixType>
class UonBasicConfiguration_ : public BasicConfiguration
{
public:
    UonBasicConfiguration_() 
    { }
    
    virtual ~UonBasicConfiguration_() 
    { }
    
    // Combines three rotational and one translation vector to one 4x4 matrix.
    virtual bool build_transformation_matrix(const VectorType& x_rotation, 
                                             const VectorType& y_rotation,
                                             const VectorType& z_rotation,
                                             const VectorType& translation,
                                             MatrixType& transformation_matrix) const = 0;

    // Parses values read by the IniReader and stores in the correspondent data fields.
    virtual bool parse();

    // Data fields
    std::vector<VectorType> contour_points_2d;
    std::vector<VectorType> contour_points_3d;
    MatrixType image_transformation;
    VectorType sub_image_in_pix_cos_ul;
    VectorType sub_image_in_pix_cos_ur;
    VectorType image_scaling;
    VectorType water_meat_compensation;

protected:

    // Stores values in the input vector by parsing a string value from an ini file.
    // A key in taken from a default section and a default subsection.
    // If no such key exists, or value parsing fails, returns false.
    virtual bool parse_vector(const String& key_name, VectorType& vector) const;

    // Stores values in the input vector by parsing a string value from an ini file.
    // A key in taken from a 'section_name' section and a 'subsection_name' subsection.
    // Key's name should match the input key_name_pattern.
    // If no such key exists, or more than one key matches the input name,
    // or value parsing fails, returns false.
    virtual bool parse_vector_by_pattern(const String& section_name,
                                         const String& subsection_name,
                                         const String& key_name_pattern,
                                         VectorType& vector) const;

    // Get the value from IniReader by it's name and trims unnecessary symbols at
    // both sides. Default section and subsection are used. If no such value found,
    // return 'def_value'.
    virtual String prepare_parsed_vector(const String& key_name,
                                         const String& def_value) const;

    // Gets the value from IniReader by it's name and trims unnecessary symbols at
    // both sides. If no such value found, returns 'def_value'.
    virtual String prepare_parsed_vector(const String& section_name,
                                         const String& subsection_name,
                                         const String& key_name,
                                         const String& def_value) const;

    // Parses a string into a vector of floats and stores it in an input vector.
    // If parsing fails, returns false and doesn't change an input vector.
    // Input string is considered to be in a standard UON vector format
    // with VECTOR_TRIM_SYMBOLS at both sides and whitespace delimiters.
    virtual bool string_to_vector(const String& str, VectorType& vector) const = 0;

}; // class UonBasicConfiguration


// Basic template class for processing UON config format.
template <class VectorType, class MatrixType>
class UonBasicConfiguration : public UonBasicConfiguration_<VectorType, MatrixType>
{ };

// UonBasicConfiguration<> class specialization for Bo's types.
template <>
class UonBasicConfiguration<Vector<float, 3>, BoFloatMatrix> :
    public UonBasicConfiguration_<Vector<float, 3>, BoFloatMatrix>
{
public:
    UonBasicConfiguration()
    {   image_transformation.resize(4, 4); 
    }

    virtual ~UonBasicConfiguration()
    { }

    virtual bool build_transformation_matrix(const Vector<float, 3>& x_rotation, 
                                             const Vector<float, 3>& y_rotation,
                                             const Vector<float, 3>& z_rotation,
                                             const Vector<float, 3>& translation,
                                             BoFloatMatrix& transformation_matrix)
                                                const;
protected:
    virtual bool string_to_vector(const String& str, Vector<float, 3>& vector) const;
}; // class UonBasicConfiguration<Vector<float, 3>, BoFloatMatrix>

#ifdef ATK_SUPPORTED
template <>
class UonBasicConfiguration<vec3<float>, mat4<float>> : 
    public UonBasicConfiguration_<vec3<float>, mat4<float>>
{
public:
    UonBasicConfiguration() 
    {   image_transformation.setIdentity();
    };
    
    virtual bool build_transformation_matrix(const vec3<float>& x_rotation, 
                                             const vec3<float>& y_rotation,
                                             const vec3<float>& z_rotation,
                                             const vec3<float>& translation,
                                             mat4<float>& transformation_matrix) const;
protected:
    virtual bool string_to_vector(const String& str, vec3<float>& vector) const;
}; // class UonBasicConfiguration<vec3<float>, mat4<float>>
#endif

// Typedefs of configurations.
typedef UonBasicConfiguration<Vector<float, 3>, BoFloatMatrix> UonBasicConfig;
#ifdef ATK_SUPPORTED
typedef UonBasicConfiguration<vec3<float>, mat4<float>> UonBasicConfigAtk;
#endif

//------------------------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------------------------
template <class VectorType, class MatrixType>
bool UonBasicConfiguration_<VectorType, MatrixType>::parse()
{
    // Read ContourPoints 2D
    contour_points_2d.clear();
    Strings contour_points_2d_names = ini_reader_.get_key_names_by_pattern(
                uon_config::CONTOUR_POINT_2D_NAME_PATTERN);
    VectorType contour_point_2d;
    for (int i = 0; i < static_cast<int>(contour_points_2d_names.size()); ++i)
    {
        if (this->parse_vector(contour_points_2d_names[i], contour_point_2d) == true)
            contour_points_2d.push_back(contour_point_2d);        
    }

    // Read ContourPoints 3D
    contour_points_3d.clear();
    Strings contour_points_3d_names = ini_reader_.get_key_names_by_pattern(
                uon_config::CONTOUR_POINT_3D_NAME_PATTERN);
    VectorType contour_point_3d;
    for (int i = 0; i < static_cast<int>(contour_points_3d_names.size()); ++i)
    {
        if (this->parse_vector(contour_points_3d_names[i], contour_point_3d) == true)
            contour_points_3d.push_back(contour_point_3d);
    }
    
    // Read ImageTransformation
    VectorType x_rotation, y_rotation, z_rotation, translation;
    this->parse_vector(uon_config::IMAGE_TRANSFORM_ROT_X_NAME, x_rotation);
    this->parse_vector(uon_config::IMAGE_TRANSFORM_ROT_Y_NAME, y_rotation);
    this->parse_vector(uon_config::IMAGE_TRANSFORM_ROT_Z_NAME, z_rotation);
    this->parse_vector(uon_config::IMAGE_TRANSFORM_TRANS_NAME, translation);
    this->build_transformation_matrix(x_rotation, y_rotation, z_rotation,
        translation, image_transformation);

    // Read SubImageInPixCosUL
    this->parse_vector(uon_config::SUB_IMAGE_IN_PIX_COS_UL_NAME, sub_image_in_pix_cos_ul);

    // Read SubImageInPixCosUR
    this->parse_vector(uon_config::SUB_IMAGE_IN_PIX_COS_UR_NAME, sub_image_in_pix_cos_ur);

    // Read ImageScaling
    this->parse_vector(uon_config::IMAGE_SCALING_NAME, image_scaling);

    // Read WaterMeatCompensation
    this->parse_vector(uon_config::WATER_MEAT_COMPENSATION_NAME, water_meat_compensation);

    return true;
}

template <class VectorType, class MatrixType>
bool UonBasicConfiguration_<VectorType, MatrixType>::
parse_vector(const String& key_name, VectorType& vector) const
{
    bool result = false;
    String default_value = "";
    String value = prepare_parsed_vector(key_name, default_value);
    if (value != default_value)
        result = this->string_to_vector(value, vector);

    return result;
}

template <class VectorType, class MatrixType>
bool UonBasicConfiguration_<VectorType, MatrixType>::
parse_vector_by_pattern(const String& section_name, const String& subsection_name, 
                        const String& key_name_pattern, VectorType& vector) const
{
    bool result = false;
    String default_value = "";
    Strings key_names = ini_reader_.get_key_names_by_pattern(section_name, 
        subsection_name, key_name_pattern);

    // Only one name should match the key_name_pattern, otherwise it is an error.
    if (key_names.size() == 1)
    {
        String value = this->prepare_parsed_vector(section_name, subsection_name,
                                                   key_names[0], default_value);
        if (value != default_value)
           result = this->string_to_vector(value, vector);
    }

    return result;
}

template <class VectorType,  class MatrixType>
String UonBasicConfiguration_<VectorType, MatrixType>::
prepare_parsed_vector(const String& key_name, const String& def_value) const
{
    String value = ini_reader_.get_value(key_name, def_value);
    boost::trim_if(value, boost::is_any_of(uon_config::VECTOR_TRIM_SYMBOLS));
    return value;
}

template <class VectorType,  class MatrixType>
String UonBasicConfiguration_<VectorType, MatrixType>::
prepare_parsed_vector(const String& section_name, const String& subsection_name,
                      const String& key_name, const String& def_value) const
{
    String value = ini_reader_.get_value(section_name, subsection_name, 
                                         key_name, def_value);
    boost::trim_if(value, boost::is_any_of(uon_config::VECTOR_TRIM_SYMBOLS));
    return value;
}

inline
bool UonBasicConfiguration<Vector<float, 3>, BoFloatMatrix>::
build_transformation_matrix(const Vector<float, 3>& x_rotation, 
                            const Vector<float, 3>& y_rotation,
                            const Vector<float, 3>& z_rotation,
                            const Vector<float, 3>& translation,
                            BoFloatMatrix& transformation_matrix) const
{
    bool result = false;

    if (x_rotation.size() == 3 && y_rotation.size() == 3 &&
        z_rotation.size() == 3 && translation.size() == 3)
    {
        transformation_matrix.resize(4, 4);

        transformation_matrix(0, 0) = x_rotation.x();
        transformation_matrix(1, 0) = x_rotation.y();
        transformation_matrix(2, 0) = x_rotation.z();
        transformation_matrix(3, 0) = 0.0;

        transformation_matrix(0, 1) = y_rotation.x();
        transformation_matrix(1, 1) = y_rotation.y();
        transformation_matrix(2, 1) = y_rotation.z();
        transformation_matrix(3, 1) = 0.0;

        transformation_matrix(0, 2) = z_rotation.x();
        transformation_matrix(1, 2) = z_rotation.y();
        transformation_matrix(2, 2) = z_rotation.z();
        transformation_matrix(3, 2) = 0.0;

        transformation_matrix(0, 3) = translation.x();
        transformation_matrix(1, 3) = translation.y();
        transformation_matrix(2, 3) = translation.z();
        transformation_matrix(3, 3) = 1.0;

        result = true;
    }

    return result;
}

inline
bool UonBasicConfiguration<Vector<float, 3>, BoFloatMatrix>::
string_to_vector(const String& str, Vector<float, 3>& vector) const
{
    float x, y, z;
    try
    {
        std::stringstream ss(str);
        ss >> x;
        ss >> y;
        ss >> z;
    }
    catch (...)
    {   return false;  
    };

    vector.x() = x;
    vector.y() = y;
    vector.z() = z;
    return true;
}

#ifdef ATK_SUPPORTED
inline
bool UonBasicConfiguration<vec3<float>, mat4<float>>::
string_to_vector(const String& str, vec3<float>& vector) const
{
    float x, y, z;
    try
    {
        std::stringstream ss(str);
        ss >> x;
        ss >> y;
        ss >> z;
    }
    catch (...)
    {   return false;  
    };

    vector.x = x;
    vector.y = y;
    vector.z = z;
    return true;
}

inline
bool UonBasicConfiguration<vec3<float>, mat4<float>>::
build_transformation_matrix(const vec3<float>& x_rotation, 
                            const vec3<float>& y_rotation,
                            const vec3<float>& z_rotation,
                            const vec3<float>& translation,
                            mat4<float>& transformation_matrix) const
{
    transformation_matrix.x = vec4<float>(x_rotation, 0);
    transformation_matrix.y = vec4<float>(y_rotation, 0);
    transformation_matrix.z = vec4<float>(z_rotation, 0);
    transformation_matrix.w = vec4<float>(translation, 1);
    return true;
}
#endif

} // namespace config_parser
} // namespace io
} // namespace bo

#endif // UON_BASIC_CONFIGURATION_HPP_7A47DDA3_838A_4177_AFD0_227E7B76FBAB_
