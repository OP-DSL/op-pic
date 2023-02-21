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

#include "fempic.h"

const int num_iterations = 50;

//*************************************************************************************************
void oppic_par_loop_inject__InjectIons(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_particle_all__MoveToCells(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_particle_inject__MoveToCells(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_all__WeightFieldsToParticles(oppic_set,oppic_arg,oppic_arg);
void oppic_par_loop_all__MoveParticles(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_all__ResetIonDensity(oppic_set,oppic_arg);
void oppic_par_loop_all__WeightParticleToMeshNodes(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    int n_nodes = 0, n_cells = 0, ts = 0;
    double dt = 1e-7, particle_remainder = 0.0;

    int *cell_to_nodes_tmp, *cell_to_cell_tmp;                                                      // cell related temporary fields arrays
    double *cell_ef_tmp,*cell_det_tmp, *cell_volume_tmp;                                            // cell related temporary fields arrays

    double *node_bnd_pot_tmp, *node_pot_tmp, *node_ion_den_tmp, *node_pos_tmp, *node_volume_tmp;    // node related temporary fields arrays

    Volume volume;
    EnrichArrays(
        volume, 
        node_bnd_pot_tmp, node_pot_tmp, node_pos_tmp, node_volume_tmp, node_ion_den_tmp,
        cell_to_nodes_tmp, cell_to_cell_tmp, cell_det_tmp, cell_volume_tmp, cell_ef_tmp
    );

    FESolver solver(volume, argc, argv);
    solver.startAssembly();
    solver.preAssembly(node_bnd_pot_tmp);

    n_nodes = volume.nodes.size();
    n_cells = volume.cells.size();

    printf("**************** n_nodes %d n_cells %d *******************\n", n_nodes, n_cells);


    {    // Start Scope for oppic
        oppic_init(argc, argv, 0);

        oppic_set nodes_set            = oppic_decl_set(n_nodes, "mesh_nodes");
        oppic_set cells_set            = oppic_decl_set(n_cells, "mesh_cells");
        oppic_set particles_set        = oppic_decl_particle_set("particles", cells_set); 

        oppic_map cell_to_nodes_map    = oppic_decl_map(cells_set, nodes_set, NODES_PER_CELL,  cell_to_nodes_tmp, "cell_to_nodes_map");   delete[] cell_to_nodes_tmp;
        oppic_map cell_to_cell_map     = oppic_decl_map(cells_set, cells_set, NEIGHBOUR_CELLS, cell_to_cell_tmp,  "cell_to_cell_map");    delete[] cell_to_cell_tmp;

        oppic_dat cell_determinants    = oppic_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), "double", sizeof(double), (char*)cell_det_tmp,    "cell_determinants");     delete[] cell_det_tmp;
        oppic_dat cell_volume          = oppic_decl_dat(cells_set, 1,                            "double", sizeof(double), (char*)cell_volume_tmp, "cell_volume");           delete[] cell_volume_tmp;
        oppic_dat cell_electric_field  = oppic_decl_dat(cells_set, DIMENSIONS,                   "double", sizeof(double), (char*)cell_ef_tmp,     "cell_electric_field");   delete[] cell_ef_tmp;

        oppic_dat node_position        = oppic_decl_dat(nodes_set, DIMENSIONS, "double", sizeof(double), (char*)node_pos_tmp,     "node_position");            delete[] node_pos_tmp;
        oppic_dat node_volume          = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_volume_tmp,  "node_volume");              delete[] node_volume_tmp;
        oppic_dat node_bnd_potential   = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_bnd_pot_tmp, "node_bnd_potential");       delete[] node_bnd_pot_tmp;
        oppic_dat node_potential       = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_pot_tmp,     "node_potential");           delete[] node_pot_tmp;
        oppic_dat node_charge_density  = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)node_ion_den_tmp, "node_charge_density");      delete[] node_ion_den_tmp;

        oppic_dat part_position        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_position");
        oppic_dat part_velocity        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_velocity");    
        oppic_dat part_electric_field  = oppic_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_electric_field");
        oppic_dat part_weights         = oppic_decl_particle_dat(particles_set, NODES_PER_CELL, "double", sizeof(double), nullptr, "part_lc");
        oppic_dat part_cell_index      = oppic_decl_particle_dat(particles_set, 1,              "int",    sizeof(int),    nullptr, "part_cell_index", true);


        auto start = std::chrono::system_clock::now();

        for (ts = 0; ts < num_iterations; ts++)
        {

        // STEP 1 - Inject Ions ************************************
            oppic_increase_particle_count(particles_set, get_num_particles_to_inject(dt, particle_remainder));
            oppic_par_loop_inject__InjectIons(
                particles_set,
                oppic_arg_dat(part_position,       OP_WRITE),
                oppic_arg_dat(part_velocity,       OP_WRITE),
                oppic_arg_dat(part_electric_field, OP_WRITE),
                oppic_arg_dat(part_weights,        OP_WRITE),
                oppic_arg_dat(part_cell_index,     OP_WRITE)
            );
            oppic_par_loop_particle_inject__MoveToCells(
                particles_set,
                oppic_arg_dat(part_position,             OP_READ),
                oppic_arg_dat(part_weights,              OP_WRITE),
                oppic_arg_dat(part_cell_index,           OP_WRITE),
                oppic_arg_dat(part_velocity,             OP_INC),
                oppic_arg_dat(cell_volume,               OP_READ, true),
                oppic_arg_dat(cell_determinants,         OP_READ, true),
                oppic_arg_dat(cell_electric_field,       OP_READ, true), 
                oppic_arg_dat(cell_to_cell_map,          OP_READ, true),
                oppic_arg_gbl(&(bool_true), 1, "bool",   OP_READ),
                oppic_arg_gbl(&dt,          1, "double", OP_READ)
            );

        // STEP X - Misc - make ion_density to zero ************************************
            oppic_par_loop_all__ResetIonDensity(
                nodes_set,
                oppic_arg_dat(node_charge_density, OP_WRITE)
            );

        // STEP 2 - Weight mesh field values to particle positions ************************************
            oppic_par_loop_all__WeightFieldsToParticles(    
                particles_set,  
                oppic_arg_dat(part_electric_field, OP_WRITE),
                oppic_arg_dat(cell_electric_field, OP_READ, true)
            );                        

        // STEP 3 - Move the particles with the influence of the fields ************************************
            oppic_par_loop_all__MoveParticles(
                particles_set,
                oppic_arg_dat(part_position,         OP_INC),
                oppic_arg_dat(part_velocity,         OP_INC),
                oppic_arg_dat(part_electric_field,   OP_READ),
                oppic_arg_gbl(&dt,  1,     "double", OP_READ)
            );
            oppic_par_loop_particle_all__MoveToCells(
                particles_set,
                oppic_arg_dat(part_position,            OP_READ),
                oppic_arg_dat(part_weights,             OP_WRITE),
                oppic_arg_dat(part_cell_index,          OP_WRITE),
                oppic_arg_dat(cell_volume,              OP_READ, true),
                oppic_arg_dat(cell_determinants,        OP_READ, true),
                oppic_arg_dat(cell_to_cell_map,         OP_READ, true),
                oppic_arg_gbl(&(bool_false), 1, "bool", OP_READ)
            );
          
        // STEP 4 - Gather the contribution from particle movement to the mesh ************************************
            oppic_par_loop_all__WeightParticleToMeshNodes(
                particles_set,
                oppic_arg_dat(part_weights,                              OP_READ),
                oppic_arg_dat(node_charge_density, 0, cell_to_nodes_map, OP_INC, true),
                oppic_arg_dat(node_charge_density, 1, cell_to_nodes_map, OP_INC, true),
                oppic_arg_dat(node_charge_density, 2, cell_to_nodes_map, OP_INC, true),
                oppic_arg_dat(node_charge_density, 3, cell_to_nodes_map, OP_INC, true),
                oppic_arg_dat(node_volume,         0, cell_to_nodes_map, OP_READ, true),
                oppic_arg_dat(node_volume,         1, cell_to_nodes_map, OP_READ, true),
                oppic_arg_dat(node_volume,         2, cell_to_nodes_map, OP_READ, true),
                oppic_arg_dat(node_volume,         3, cell_to_nodes_map, OP_READ, true)
            );

        // STEP 5 - Solve field values on the mesh points ************************************
            //matrix | TODO : Try to make these in to kernels and use oppic_par_loop
            solver.SolveFields(
                node_charge_density,    // OP_READ
                node_potential,         // OP_WRITE
                node_bnd_potential,     // OP_READ
                cell_electric_field     // OP_WRITE
            );

        // STEP 6 - Log and/or print values to files ************************************ 
            if (print_all_to_file || ((ts + 1) == num_iterations))
            {
                std::string f = std::string("F_") + std::to_string(ts + 1);

                oppic_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");
                oppic_print_dat_to_txtfile(node_potential, f.c_str(), "node_potential.dat");
                oppic_print_dat_to_txtfile(cell_electric_field, f.c_str(), "cell_electric_field.dat");

                oppic_print_dat_to_txtfile(part_position, f.c_str(), "part_position.dat");
                oppic_print_dat_to_txtfile(part_velocity, f.c_str(), "part_velocity.dat");
                oppic_print_dat_to_txtfile(part_weights, f.c_str(), "part_weights.dat");
                oppic_print_dat_to_txtfile(part_cell_index, f.c_str(), "part_cell_index.dat");
            }
            std::cout << "ts: " << ts << "\t num_particles: " << particles_set->size << std::endl;

        } // end of time marching loop


        std::chrono::duration<double> diff = std::chrono::system_clock::now() - start;
        std::cout << "\nFEMPIC - Time to iterate " << ts << " takes <chrono>: " << diff.count() << " s\n\n";

        oppic_exit();

    } // End Scope for oppic


    return 0;
}

//*************************************************************************************************



