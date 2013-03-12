
/******************************************************************************

  I/O for Mesh class. RPly library is used for working with .ply files.

  Copyright (c) 2010 - 2012
  Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef MESH_IO_HPP_C536C1D7_8F7D_4396_90C3_84DFCB3902C5_
#define MESH_IO_HPP_C536C1D7_8F7D_4396_90C3_84DFCB3902C5_

#include <string>

#include "bo/config.hpp"
#include "bo/mesh.hpp"

namespace bo {
namespace io {

// IO functions, allow to read mesh from and write to a .ply files. Currently, only
// meshes with float vertices are supported, moreover no checks are performed to
// determine, whether the .ply file has actually floats. That behaviour should be
// changed. At least two different formats should be supported: integral and double.
bo::Mesh<float> BO_DECL mesh_from_ply(const std::string& file_path);
bool BO_DECL mesh_to_ply(const bo::Mesh<float>& mesh, const std::string& file_path);

} // namespace io
} // namespace bo

#endif // MESH_IO_HPP_C536C1D7_8F7D_4396_90C3_84DFCB3902C5_
