
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

#include "fempic.h"
#include "../lib_oppic/oppic.h"
#include "kernels.h"

const int num_iterations = 50;

//********************************************* MAIN ****************************************************
int main() 
{
    int n_nodes = 0, n_cells = 0, ts = 0;
    double dt = 1e-7, particle_remainder = 0.0;

    int *cell_to_nodes_tmp, *cell_to_cell_tmp;                                                      // cell related temporary fields arrays
    double *cell_ef_tmp,*cell_det_tmp, *cell_volume_tmp;                                            // cell related temporary fields arrays

    double *node_bnd_pot_tmp, *node_pot_tmp, *node_ion_den_tmp, *node_pos_tmp, *node_volume_tmp;     // node related temporary fields arrays

    Volume volume;
    EnrichArrays(volume, 
        node_bnd_pot_tmp, node_pot_tmp, node_pos_tmp, node_volume_tmp, node_ion_den_tmp,
        cell_to_nodes_tmp, cell_to_cell_tmp, cell_det_tmp, cell_volume_tmp, cell_ef_tmp
    );

    FESolver solver(volume);
    solver.startAssembly();
    solver.preAssembly(node_bnd_pot_tmp);
    
    n_nodes = volume.nodes.size();
    n_cells = volume.cells.size();

    printf("**************** n_nodes %d n_cells %d *******************\n", n_nodes, n_cells);


    {    // Start Scope for oppic
        op_init();

        op_set nodes_set            = op_decl_set(n_nodes, "mesh_nodes");
        op_set cells_set            = op_decl_set(n_cells, "mesh_cells");
        op_set particles_set        = op_decl_particle_set("particles", cells_set); 

        op_map cell_to_nodes_map    = op_decl_map(cells_set, nodes_set, NODES_PER_CELL,  cell_to_nodes_tmp, "cell_to_nodes_map"); delete[] cell_to_nodes_tmp;
        op_map cell_to_cell_map     = op_decl_map(cells_set, cells_set, NEIGHBOUR_CELLS, cell_to_cell_tmp,  "cell_to_cell_map");  delete[] cell_to_cell_tmp;

        op_dat cell_determinants    = op_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), "double", sizeof(double), (char*)cell_det_tmp,    "cell_determinants");   delete[] cell_det_tmp;
        op_dat cell_volume          = op_decl_dat(cells_set, 1,                            "double", sizeof(double), (char*)cell_volume_tmp, "cell_volume");         delete[] cell_volume_tmp;
        op_dat cell_electric_field  = op_decl_dat(cells_set, DIMENSIONS,                   "double", sizeof(double), (char*)cell_ef_tmp,     "cell_electric_field"); delete[] cell_ef_tmp;

        op_dat node_position        = op_decl_dat(nodes_set, DIMENSIONS, "double", sizeof(double), (char*)node_pos_tmp,        "node_position");       delete[] node_pos_tmp;
        op_dat node_volume          = op_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_volume_tmp,     "node_volume");         delete[] node_volume_tmp;
        op_dat node_bnd_potential   = op_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_bnd_pot_tmp,    "node_bnd_potential");  delete[] node_bnd_pot_tmp;
        op_dat node_potential       = op_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_pot_tmp,        "node_potential");      delete[] node_pot_tmp;
        op_dat node_charge_density  = op_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_ion_den_tmp,    "node_charge_density"); delete[] node_ion_den_tmp;

        op_dat part_position        = op_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_position");
        op_dat part_velocity        = op_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_velocity");    
        op_dat part_electric_field  = op_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_electric_field");
        op_dat part_weights         = op_decl_particle_dat(particles_set, NODES_PER_CELL, "double", sizeof(double), nullptr, "part_lc");
        op_dat part_cell_index      = op_decl_particle_dat(particles_set, 1,              "int",    sizeof(int),    nullptr, "part_cell_index", true);


        auto start = std::chrono::system_clock::now();

        for (ts = 0; ts < num_iterations; ts++)
        {

        // STEP 1 - Inject Ions ************************************
            op_increase_particle_count(particles_set, get_num_particles_to_inject(dt, particle_remainder));
            op_par_loop( 
                inject_ions__kernel, "inject_ions__kernel",
                particles_set, OP_ITERATE_INJECTED,
                op_arg_dat(part_position,       OP_WRITE),
                op_arg_dat(part_velocity,       OP_WRITE),
                op_arg_dat(part_electric_field, OP_WRITE),
                op_arg_dat(part_weights,        OP_WRITE),
                op_arg_dat(part_cell_index,     OP_WRITE)
            );
            op_par_loop_particle(
                move_particle_to_cell__kernel, "move_particle_to_cell__kernel",
                particles_set, OP_ITERATE_INJECTED,
                op_arg_dat(part_position,           OP_WRITE),
                op_arg_dat(part_weights,            OP_WRITE),
                op_arg_dat(part_cell_index,         OP_WRITE),
                op_arg_dat(cell_volume,             OP_READ, true),
                op_arg_dat(cell_determinants,       OP_READ, true),
                op_arg_dat(cell_to_cell_map,        OP_READ, true),
                op_arg_gbl(&(bool_true), 1, "bool", OP_READ)
            );
            op_par_loop(
                enrich_velocity__kernel, "enrich_velocity__kernel",
                particles_set, OP_ITERATE_INJECTED,
                op_arg_dat(part_velocity,       OP_WRITE),
                op_arg_dat(cell_electric_field, OP_READ, true),
                op_arg_gbl(&dt, 1, "double",    OP_READ)
            );
            op_particle_sort(particles_set);

        // STEP X - Misc - make ion_density to zero ************************************
            op_par_loop(
                reset_ion_density__kernel, "reset_ion_density__kernel",
                nodes_set, OP_ITERATE_ALL,
                op_arg_dat(node_charge_density, OP_WRITE)
            );

        // STEP 2 - Weight mesh field values to particle positions ************************************
            op_par_loop(
                weight_fields_to_particles__kernel, "weight_fields_to_particles__kernel",    
                particles_set, OP_ITERATE_ALL,
                op_arg_dat(part_electric_field, OP_WRITE),
                op_arg_dat(cell_electric_field, OP_READ, true)
            );

        // STEP 3 - Move the particles with the influence of the fields ************************************
            op_par_loop(
                move_particles__kernel, "move_particles__kernel",
                particles_set, OP_ITERATE_ALL,
                op_arg_dat(part_position,         OP_INC),
                op_arg_dat(part_velocity,         OP_INC),
                op_arg_dat(part_electric_field,   OP_READ),
                op_arg_gbl(&dt,  1,     "double", OP_READ)
            );
            op_par_loop_particle(
                move_particle_to_cell__kernel, "move_particle_to_cell__kernel",
                particles_set, OP_ITERATE_ALL,
                op_arg_dat(part_position,            OP_WRITE),
                op_arg_dat(part_weights,             OP_WRITE),
                op_arg_dat(part_cell_index,          OP_WRITE),
                op_arg_dat(cell_volume,              OP_READ, true),
                op_arg_dat(cell_determinants,        OP_READ, true),
                op_arg_dat(cell_to_cell_map,         OP_READ, true),
                op_arg_gbl(&(bool_false), 1, "bool", OP_READ)
            );
            op_particle_sort(particles_set);

        // STEP 4 - Gather the contribution from particle movement to the mesh ************************************
            op_par_loop(
                weight_particle_to_mesh_nodes__kernel, "weight_particle_to_mesh_nodes__kernel",
                particles_set, OP_ITERATE_ALL,
                op_arg_dat(part_weights,                              OP_READ),
                op_arg_dat(part_cell_index,                           OP_READ),
                op_arg_dat(node_charge_density, 0, cell_to_nodes_map, OP_INC, true),
                op_arg_dat(node_charge_density, 1, cell_to_nodes_map, OP_INC, true),
                op_arg_dat(node_charge_density, 2, cell_to_nodes_map, OP_INC, true),
                op_arg_dat(node_charge_density, 3, cell_to_nodes_map, OP_INC, true),
                op_arg_dat(node_volume,         0, cell_to_nodes_map, OP_READ, true),
                op_arg_dat(node_volume,         1, cell_to_nodes_map, OP_READ, true),
                op_arg_dat(node_volume,         2, cell_to_nodes_map, OP_READ, true),
                op_arg_dat(node_volume,         3, cell_to_nodes_map, OP_READ, true)
            );

        // STEP 5 - Solve field values on the mesh points ************************************
            //matrix | TODO : Try to make these in to kernels and use op_par_loop 
            solver.SolveFields(                            /*field solve*/
                (double *)node_charge_density->data,
                (double *)node_potential->data,
                (double *)node_bnd_potential->data, 
                (double *)cell_electric_field->data            
            );    

        // STEP 6 - Log and/or print values to files ************************************
            if (print_all_to_file || ((ts + 1) == num_iterations))
            {
                std::string f = std::string("F_") + std::to_string(ts + 1);
                print_fields_m((double*)(node_charge_density->data), (double*)(node_potential->data), (double*)(cell_electric_field->data), n_nodes, n_cells, f.c_str());
                print_particles_m((double*)(part_position->data), (double*)(part_velocity->data), (int*)(part_cell_index->data), particles_set->size, f.c_str());
            }
            std::cout << "ts: " << ts << "\t num_particles: " << particles_set->size << std::endl;

        } // end of time marching loop


        std::chrono::duration<double> diff = std::chrono::system_clock::now() - start;
        std::cout << "\nFEMPIC - Time to iterate " << ts << " takes <chrono>: " << diff.count() << " s\n\n";

        op_exit();

    } // End Scope for oppic


    return 0;
}

//*************************************************************************************************