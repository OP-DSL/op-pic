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

// USER WRITTEN CODE

#include "simpic.h"
#include "../lib_oppic/oppic.h"
#include "kernels.h"

const bool print_all_to_file	= false;

//*********************************************MAIN****************************************************
int main( int argc, char* argv[] ) 
{
	parsecmdline(argc, argv);
	init();

	int n_cells = nc, n_nodes = ng, n_particles  = npart, loop_count = 0;
	double max_time = t, time_step = dt;

	int *cell_to_nodes_tmp, *node_index_tmp, *part_cell_index_tmp;
	double *node_field_E_tmp, *node_field_J_tmp, *node_field_P_tmp, *node_xlocal_tmp, *part_position_x_tmp, *part_velocity_x_tmp, *part_field_E_tmp;

	init_fields(node_field_E_tmp, node_field_J_tmp, node_field_P_tmp, node_xlocal_tmp, node_index_tmp, cell_to_nodes_tmp);
	init_particles(part_position_x_tmp, part_velocity_x_tmp, part_field_E_tmp, part_cell_index_tmp);

	printf("**************** n_nodes %d n_cells %d *******************\n", n_nodes, n_cells);

		
	{	// Start Scope for oppic
		op_init();

		op_set nodes_set  			= op_decl_set(n_nodes, "mesh_nodes");
		op_set cells_set  			= op_decl_set(n_cells, "mesh_cells");
		op_set particles_set   		= op_decl_particle_set(n_particles, "particles", cells_set);  
	
		op_map cell_to_nodes_map 	= op_decl_map(cells_set, nodes_set, 2, cell_to_nodes_tmp, "cell_to_edges_map");  delete[] cell_to_nodes_tmp; // LHS node, RHS cell
		
		op_dat node_field_E 		= op_decl_dat(nodes_set, 1, "double",  sizeof(double),  (char*)node_field_E_tmp,  "node_field_E");	delete[] node_field_E_tmp;
		op_dat node_field_J 		= op_decl_dat(nodes_set, 1, "double",  sizeof(double),  (char*)node_field_J_tmp,  "node_field_J");	delete[] node_field_J_tmp;
		op_dat node_field_P 		= op_decl_dat(nodes_set, 1, "double",  sizeof(double),  (char*)node_field_P_tmp,  "node_field_P");	delete[] node_field_P_tmp;
		op_dat node_xlocal  		= op_decl_dat(nodes_set, 1, "double",  sizeof(double),  (char*)node_xlocal_tmp,   "node_xlocal");	delete[] node_xlocal_tmp;

		op_dat part_position_x  	= op_decl_particle_dat(particles_set,  1, "double", sizeof(double), (char*)part_position_x_tmp, "part_position_x");			delete[] part_position_x_tmp;
		op_dat part_velocity_x  	= op_decl_particle_dat(particles_set,  1, "double", sizeof(double), (char*)part_velocity_x_tmp, "part_velocity_x");			delete[] part_velocity_x_tmp;
		op_dat part_field_E	 		= op_decl_particle_dat(particles_set,  1, "double", sizeof(double), (char*)part_field_E_tmp,	"part_field_E");		 	delete[] part_field_E_tmp;
		op_dat part_cell_index  	= op_decl_particle_dat(particles_set,  1, "int",	sizeof(int),	(char*)part_cell_index_tmp, "part_cell_index", true);	delete[] part_cell_index_tmp;


		auto start = std::chrono::system_clock::now();

		for (double tt = 0; tt < max_time; tt += time_step, ++loop_count)
		{

		// STEP X - Misc - make current_density to zero ************************************
			op_par_loop(
				reset_current_density__kernel, "ResetCurrentDensity",
				nodes_set, OP_ITERATE_ALL,
				op_arg_dat(node_field_J, OP_WRITE)
			);

		// STEP 2 - Weight mesh field values to particle positions ************************************
			op_par_loop(
				weight_fields_to_particles__kernel, "WeightFieldsToParticles",
				particles_set, OP_ITERATE_ALL,
				op_arg_dat(node_field_E, 0, cell_to_nodes_map, 	OP_READ, true),  // lhs node_field_E
				op_arg_dat(node_field_E, 1, cell_to_nodes_map, 	OP_READ, true),  // rhs node_field_E  
				op_arg_dat(part_position_x, 				OP_READ),	
				op_arg_dat(part_field_E,					OP_WRITE)	 
			);

		// STEP 3 - Move the particles with the influence of the fields ************************************
			op_par_loop_particle(
				particle_push__kernel, "PushParticles",
				particles_set, OP_ITERATE_ALL,
				op_arg_dat(part_field_E,	OP_READ),
				op_arg_dat(part_velocity_x, OP_RW),
				op_arg_dat(part_position_x, OP_RW),
				op_arg_dat(part_cell_index, OP_WRITE)
			);
			//op_particle_sort(particles_set);

		// STEP 4 - Gather the contribution from particle movement to the mesh ************************************
			op_par_loop(
				weight_particles_to_fields__kernel, "WeightParticlesToFields",
				particles_set, OP_ITERATE_ALL,
				op_arg_dat(node_field_J, 0, cell_to_nodes_map,	OP_INC),  // lhs node_field_J
				op_arg_dat(node_field_J, 1, cell_to_nodes_map,	OP_INC),  // rhs node_field_J
				op_arg_dat(part_position_x, 					OP_READ)  
			);
		
		// STEP 5 - Solve field values on the mesh points ************************************
			seq_field_solve_poissons_equation( // TODO : This needs to be converted to kernels
				nodes_set->size,
				(double*)node_field_J->data,
				(double*)node_field_P->data
			);

			op_par_loop(
				field_solve_sum_laplace__kernel, "FieldSolveSumLaplace",
				nodes_set, OP_ITERATE_ALL,
				op_arg_dat(node_xlocal,  OP_READ),
				op_arg_dat(node_field_P, OP_WRITE)
			);

			seq_field_solve_get_potential_gradient( // TODO : This needs to be converted to kernels
				nodes_set->size,
				(double*)node_field_E->data,
				(double*)node_field_P->data
			);

		// STEP 6 - Log and/or print values to files ************************************ 
			if (print_all_to_file || ((tt + time_step) > max_time))
			{   
				std::string f = std::string("F_") + std::to_string(loop_count);
				printDataToFiles(n_nodes, n_particles, 
								(double*)node_field_E->data, (double*)node_field_J->data, (double*)node_field_P->data, 
								(double*)part_position_x->data, (double*)part_velocity_x->data, (int*)part_cell_index->data, f);
			} 

			printf("Looping dt = %+2.15lE, tmax = %+2.15lE loop %d ******************\n", tt, max_time, loop_count); 
		}

		std::chrono::duration<double> diff = std::chrono::system_clock::now() - start;
		std::cout << "\nFEMPIC - Time to iterate " << loop_count << " takes <chrono>: " << diff.count() << " s\n\n";

		op_exit();

	} // End Scope for oppic


	return 0;
}

