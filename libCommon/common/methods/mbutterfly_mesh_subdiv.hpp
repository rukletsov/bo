
/******************************************************************************

  mbutterfly_mesh_subdiv.hpp, v 1.0.0 2012.02.09

  Implementation of the Modified Butterfly subdivision surface method.

  Reference paper: Zorin, D., Schroeder, P. and Sweldens, W., 
  “Interpolating subdivision for meshes with arbitrary topology”, 
  Computer Graphics, Ann. Conf. Series 30, 189-192 (1996).

  Copyright (c) 2012
  Dzmitry Hlindzich <hlindzich@gmail.com>
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

#ifndef MBUTTERFLY_MESH_SUBDIV_HPP_0CE908A8_0799_4C4E_92EC_3DCDAB62386E_
#define MBUTTERFLY_MESH_SUBDIV_HPP_0CE908A8_0799_4C4E_92EC_3DCDAB62386E_

#include <common/mesh.hpp>

namespace common {
namespace methods {

namespace surfaces {

    // Performs the given number of subdivision iterations on the source mesh
    // using the Modified Butterfly surface subdivision scheme (see the reference 
    // paper) and returns the resulted mesh.
    Mesh mbutterfly_subdivision(const Mesh &source, int iterations);

} //namespace surfaces

} //namespace methods
} //namespace common

#endif //MBUTTERFLY_MESH_SUBDIV_HPP_0CE908A8_0799_4C4E_92EC_3DCDAB62386E_
