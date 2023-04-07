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

//*********************************************
// AUTO GENERATED CODE
//*********************************************


#include "../fempic.h"
#include <oppic_hip.h>

#define GPU_THREADS_PER_BLOCK 128

// TODO : This should be removed
double CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0;

//****************************************
__constant__ double CONST_spwt_cuda, CONST_ion_velocity_cuda = 0, CONST_dt_cuda = 0, CONST_plasma_den_cuda = 0, CONST_mass_cuda = 0, CONST_charge_cuda = 0;
void oppic_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_spwt"))              cutilSafeCall(hipMemcpyToSymbol(CONST_spwt_cuda, data, dim*size));
    else if (!strcmp(name,"CONST_ion_velocity")) cutilSafeCall(hipMemcpyToSymbol(CONST_ion_velocity_cuda, data, dim*size));
    else if (!strcmp(name,"CONST_dt"))           cutilSafeCall(hipMemcpyToSymbol(CONST_dt_cuda, data, dim*size));
    else if (!strcmp(name,"CONST_plasma_den"))   cutilSafeCall(hipMemcpyToSymbol(CONST_plasma_den_cuda, data, dim*size));
    else if (!strcmp(name,"CONST_mass"))         cutilSafeCall(hipMemcpyToSymbol(CONST_mass_cuda, data, dim*size));
    else if (!strcmp(name,"CONST_charge"))       cutilSafeCall(hipMemcpyToSymbol(CONST_charge_cuda, data, dim*size));
    else std::cerr << "error: unknown const name" << std::endl;

    // TODO : This block should be removed
    {    if (!strcmp(name,"CONST_spwt"))              CONST_spwt = *((double*)data);
        else if (!strcmp(name,"CONST_ion_velocity")) CONST_ion_velocity = *((double*)data);
        else if (!strcmp(name,"CONST_dt"))           CONST_dt = *((double*)data);
        else if (!strcmp(name,"CONST_plasma_den"))   CONST_plasma_den = *((double*)data);
        else if (!strcmp(name,"CONST_mass"))         CONST_mass = *((double*)data);
        else if (!strcmp(name,"CONST_charge"))       CONST_charge = *((double*)data);
        else std::cerr << "error: unknown const name" << std::endl; } // TODO : This block should be removed
}
//****************************************

//*************************************************************************************************
#include "oppic_inject__Increase_particle_count.cpp"

//*************************************************************************************************
#include "oppic_par_loop_inject__InjectIons.cpp"

//*************************************************************************************************
#include "oppic_par_loop_particle_all__MoveToCells.cpp"

//*************************************************************************************************
#include "oppic_par_loop_all__ComputeNodeChargeDensity.cpp"

//*************************************************************************************************
#include "oppic_par_loop_all__ComputeElectricField.cpp"

//*************************************************************************************************