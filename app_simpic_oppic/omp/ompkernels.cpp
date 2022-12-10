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

#include "../../lib_oppic/oppic.h"
#include "../kernels.h"

//*************************************************************************************************
void op_par_loop_all__ResetCurrentDensity(
	op_set set, 	// nodes_set
	op_arg arg0		// node_charge_density
	)
{ TRACE_ME;
	
	if (OP_DEBUG) printf("SIMPIC - op_par_loop_all__ResetCurrentDensity [%d]\n", set->size);
	
}

//*************************************************************************************************
void op_par_loop_all__WeightFieldsToParticles(
	op_set set, 	// particles
	op_arg arg0,	// lhs node0_field_E
	op_arg arg1,	// rhs node1_field_E
	op_arg arg2,	// particle0_position_x
	op_arg arg3  	// particle0_field_E
	)
{ TRACE_ME;

	if (OP_DEBUG) printf("SIMPIC - op_par_loop_all__WeightFieldsToParticles [%d]\n", set->size);

}

//*************************************************************************************************
void op_par_loop_particle_all__PushParticles( 
	op_set set, 	// particles
	op_arg arg0,	// particle0_field_E
	op_arg arg1,	// particle0_velocity_x
	op_arg arg2,	// particle0_position_x
	op_arg arg3		// particle0_cell_index
	)
{ TRACE_ME;

	if (OP_DEBUG) printf("SIMPIC - op_par_loop_particle_all__PushParticles [%d]\n", set->size);

}

//*************************************************************************************************
void op_par_loop_all__WeightParticlesToFields(
	op_set set, 	// particles
	op_arg arg0,	// lhs node0_field_J
	op_arg arg1,	// rhs node1_field_J
	op_arg arg2		// particle0_position_x
	)
{ TRACE_ME;

	if (OP_DEBUG) printf("SIMPIC - op_par_loop_all__WeightParticlesToFields [%d]\n", set->size);

}

//*************************************************************************************************
void op_par_loop_all__FieldSolveSumLaplace(
	op_set set,   // nodes
	op_arg arg0,  // node0_xlocal
	op_arg arg1)  // node0_field_P
{ TRACE_ME;

	if (OP_DEBUG) printf("SIMPIC - op_par_loop__FieldSolveSumLaplace [%d]\n", set->size);

}

//*************************************************************************************************