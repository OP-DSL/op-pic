/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// AUTO GENERATED CODE

#include "../fempic.h"
#include <oppic_cuda.h>

#define GPU_THREADS_PER_BLOCK 16

__constant__ double OP_CONST_CUDA_charge = 1.602e-19;                // TODO : Make this OP2 constants
__constant__ double OP_CONST_CUDA_mass   = (16 * 1.660538921e-27);   // TODO : Make this OP2 constants
__constant__ double OP_CONST_CUDA_spwt   = 2e2;                      // TODO : Make this OP2 constants

//*************************************************************************************************

#include "oppic_par_loop__InjectIons.cu"

#include "oppic_par_loop_particle_inject__MoveToCells.cu"

#include "oppic_par_loop__WeightFieldsToParticles.cu"

#include "oppic_par_loop__MoveParticles.cu"

#include "oppic_par_loop_particle_all__MoveToCells.cu"

#include "oppic_par_loop__ResetIonDensity.cu"

#include "oppic_par_loop__WeightParticleToMeshNodes.cu"

//*************************************************************************************************