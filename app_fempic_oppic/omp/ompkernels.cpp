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

#include "oppic_omp.h"
#include "../kernels.h"


//*************************************************************************************************
void oppic_par_loop_inject__InjectIons(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_velocity,
    oppic_arg arg2,     // part_electric_field,
    oppic_arg arg3,     // part_weights,
    oppic_arg arg4      // part_cell_index,
    )
{ TRACE_ME;
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_inject__InjectIons num_particles %d diff %d\n", set->size, set->diff);

    int nthreads = omp_get_max_threads();

    if (set->size > 0) 
    {
        //#pragma omp parallel for // CANNOT USE PARALLEL FOR NOW: Since random generator will give unexpected results
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = (set->size - set->diff) + (set->diff* thr)/nthreads;
            int finish = (set->size - set->diff) + (set->diff*(thr+1))/nthreads;
            
            for (int n = start; n < finish; n++)
            {    
                inject_ions__kernel(    
                    &((double *)arg0.data)[n * arg0.dim],    // part_position,
                    &((double *)arg1.data)[n * arg1.dim],    // part_velocity,
                    &((double *)arg2.data)[n * arg2.dim],    // part_electric_field,
                    &((double *)arg3.data)[n * arg3.dim],    // part_weights,
                    &((int *)arg4.data)[n * arg4.dim]        // part_cell_index,
                );
            }
        }
    }
}

//*************************************************************************************************
void oppic_par_loop_particle_inject__MoveToCells(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_weights,
    oppic_arg arg2,     // part_cell_index,
    oppic_arg arg3,     // part_velocity,
    oppic_arg arg4,     // cell_volume,
    oppic_arg arg5,     // cell_det,
    oppic_arg arg6,     // cell_electric_field,
    oppic_arg arg7,     // cell_connectivity_map,
    oppic_arg arg8,     // particles_injected,
    oppic_arg arg9      // dt,
    )
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_particle_inject__MoveToCells num_particles %d diff %d\n", set->size, set->diff);

    const int num_cells    = set->cells_set->size; 
    int nthreads = omp_get_max_threads();
    
    oppic_init_particle_move(set);

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        int start  = (set->size - set->diff) + (set->diff * thr)/nthreads;
        int finish = (set->size - set->diff) + (set->diff * (thr+1))/nthreads;
        
        for (int n = start; n < finish; n++)
        {        
            int& map0idx    = ((int *)set->cell_index_dat->data)[n * set->cell_index_dat->dim];    // this is the cell_index
            int move_status = (int)NEED_MOVE;
            
            do
            {
                move_injected_particles_to_cell__kernel(
                    &(move_status),
                    &((double *)arg0.data)[n * arg0.dim],        // part_pos
                    &((double *)arg1.data)[n * arg1.dim],        // part_weights
                    &((int *)arg2.data)[n * arg2.dim],           // part_cell_index
                    &((double *)arg3.data)[n * arg3.dim],        // part_velocity
                    &((double*)arg4.data)[map0idx * arg4.dim],   // cell_volume
                    &((double*)arg5.data)[map0idx * arg5.dim],   // cell_det
                    &((double*)arg6.data)[map0idx * arg6.dim],   // cell_electric_field
                    &((int*)arg7.data)[map0idx * arg7.dim],      // cell_connectivity
                    (bool *)arg8.data,                           // full_mesh_search
                    (double *)arg9.data                          // dt   
                ); 
                
            } while ((move_status == (int)NEED_MOVE) && (map0idx < num_cells));

            oppic_mark_particle_to_move(set, n, move_status);
        }
    }

    oppic_finalize_particle_move(set);
}

//*************************************************************************************************
void oppic_par_loop_all__ResetIonDensity(
    oppic_set set,     // nodes_set
    oppic_arg arg0     // node_charge_density
    )
{ TRACE_ME;
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__ResetIonDensity num_nodes %d\n", set->size);

    int nthreads = omp_get_max_threads();

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = (set->size* thr)/nthreads;
            int finish = (set->size*(thr+1))/nthreads;

            for (int n = start; n < finish; n++)
            {
                reset_ion_density__kernel(
                    &((double*)arg0.data)[n * arg0.dim]
                );
            }
        }
    }    
}

//*************************************************************************************************
void oppic_par_loop_all__WeightFieldsToParticles(
    oppic_set set,             // particles_set
    oppic_arg arg0,            // particle_ef
    oppic_arg arg1             // cell_electric_field
    )
{ TRACE_ME;
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__WeightFieldsToParticles num_particles %d\n", set->size);

    int nthreads = omp_get_max_threads();

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = (set->size* thr)/nthreads;
            int finish = (set->size*(thr+1))/nthreads;
            
            for (int n = start; n < finish; n++)
            {
                int map0idx    = ((int *)set->cell_index_dat->data)[n * set->cell_index_dat->dim];
                
                weight_fields_to_particles__kernel(
                    &((double*)arg0.data)[n * arg0.dim],            // particle_ef
                    &((double*)arg1.data)[map0idx * arg1.dim]       // cell_electric_field
                );
            }
        }
    }
}

//*************************************************************************************************
/*updates ion velocities and positions*/
void oppic_par_loop_all__MoveParticles(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_velocity,
    oppic_arg arg2,     // part_electric_field,
    oppic_arg arg3      // const dt 
    )
{ TRACE_ME;    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__MoveParticles num_particles %d\n", set->size);

    int nthreads = omp_get_max_threads();

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = (set->size* thr)/nthreads;
            int finish = (set->size*(thr+1))/nthreads;
            
            for (int n = start; n < finish; n++)
            {
                move_particles__kernel(
                    &((double *)arg0.data)[n * arg0.dim],            // part_position,
                    &((double *)arg1.data)[n * arg1.dim],            // part_velocity,
                    &((double *)arg2.data)[n * arg2.dim],            // part_electric_field,
                    (double*)arg3.data                               // const dt 
                );
            }
        }
    }
}

//*************************************************************************************************
void oppic_par_loop_particle_all__MoveToCells(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_weights,
    oppic_arg arg2,     // part_cell_index,
    oppic_arg arg3,     // cell_volume,
    oppic_arg arg4,     // cell_det,
    oppic_arg arg5,     // cell_connectivity_map,
    oppic_arg arg6      // particles_injected
    )
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_particle_all__MoveToCells num_particles %d diff %d\n", set->size, set->diff);

    const int num_cells    = set->cells_set->size; 
    int nthreads = omp_get_max_threads();
    
    oppic_init_particle_move(set);

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        int start  = (set->size * thr)/nthreads;
        int finish = (set->size * (thr+1))/nthreads;
        
        for (int n = start; n < finish; n++)
        {        
            int& map0idx    = ((int *)set->cell_index_dat->data)[n * set->cell_index_dat->dim];    // this is the cell_index
            int move_status = (int)NEED_MOVE;

            do
            {
                move_all_particles_to_cell__kernel(
                    &(move_status),
                    &((double *)arg0.data)[n * arg0.dim],         // part_pos
                    &((double *)arg1.data)[n * arg1.dim],         // part_weights
                    &((int *)arg2.data)[n * arg2.dim],            // part_cell_index
                    &((double*)arg3.data)[map0idx * arg3.dim],    // cell_volume
                    &((double*)arg4.data)[map0idx * arg4.dim],    // cell_det
                    &((int*)arg5.data)[map0idx * arg5.dim],       // cell_connectivity
                    (bool *)arg6.data                             // full_mesh_search     
                );                                                    
                
            } while ((move_status == (int)NEED_MOVE) && (map0idx < num_cells));

            oppic_mark_particle_to_move(set, n, move_status);
        }
    }

    oppic_finalize_particle_move(set);
}

//*************************************************************************************************
void oppic_par_loop_all__WeightParticleToMeshNodes(
    oppic_set set,         // particles_set
    oppic_arg arg0,        // particle_lc
    oppic_arg arg1,        // node_charge_density
    oppic_arg arg2,        // node_charge_density
    oppic_arg arg3,        // node_charge_density
    oppic_arg arg4,        // node_charge_density
    oppic_arg arg5,        // node_volumes
    oppic_arg arg6,        // node_volumes
    oppic_arg arg7,        // node_volumes
    oppic_arg arg8         // node_volumes        
    )
{ TRACE_ME;
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__WeightParticleToMeshNodes num_particles %d\n", set->size);

    int nthreads = omp_get_max_threads();

    oppic_create_thread_level_data<double>(arg1, 0.0); // For now, this does not work for particle dats due to resizing

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = (set->size* thr)/nthreads;
            int finish = (set->size*(thr+1))/nthreads;

            for (int n = start; n < finish; n++)
            {
                int map0idx    = ((int *)set->cell_index_dat->data)[n * set->cell_index_dat->dim];

                const int map1idx = arg1.map_data[map0idx * arg1.dim + 0];
                const int map2idx = arg1.map_data[map0idx * arg1.dim + 1];
                const int map3idx = arg1.map_data[map0idx * arg1.dim + 2];
                const int map4idx = arg1.map_data[map0idx * arg1.dim + 3];

                const int map5idx = arg5.map_data[map0idx * arg5.dim + 0];
                const int map6idx = arg5.map_data[map0idx * arg5.dim + 1];
                const int map7idx = arg5.map_data[map0idx * arg5.dim + 2];
                const int map8idx = arg5.map_data[map0idx * arg5.dim + 3];

                weight_particle_to_mesh_nodes__kernel(
                    &((double*)arg0.data)[n * arg0.dim],                    // part_lc
                    &((double*)arg2.dat->thread_data->at(thr))[map1idx],    // node_charge_den0
                    &((double*)arg2.dat->thread_data->at(thr))[map2idx],    // node_charge_den1
                    &((double*)arg2.dat->thread_data->at(thr))[map3idx],    // node_charge_den2
                    &((double*)arg2.dat->thread_data->at(thr))[map4idx],    // node_charge_den3
                    &((double*)arg5.data)[map5idx],                         // node_volume0
                    &((double*)arg5.data)[map6idx],                         // node_volume1
                    &((double*)arg5.data)[map7idx],                         // node_volume2
                    &((double*)arg5.data)[map8idx]                          // node_volume3
                );
            }
        }
    }

    oppic_reduce_thread_level_data<double>(arg1); // For now, this does not work for particle dats due to resizing
}

// GIVES CORRECT RESULTS WITH BELOW, THE ABOVE HAS A PROBLEM WITH NODE_CHARGE_DENSITY
// //*************************************************************************************************
// void oppic_par_loop_all__WeightParticleToMeshNodes(
//     oppic_set set,         // particles_set
//     oppic_arg arg0,        // particle_lc
//     oppic_arg arg1,        // node_charge_density
//     oppic_arg arg2,        // node_charge_density
//     oppic_arg arg3,        // node_charge_density
//     oppic_arg arg4,        // node_charge_density
//     oppic_arg arg5,        // node_volumes
//     oppic_arg arg6,        // node_volumes
//     oppic_arg arg7,        // node_volumes
//     oppic_arg arg8         // node_volumes        
//     )
// { TRACE_ME;

//     if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__WeightParticleToMeshNodes num_particles %d\n", set->size);

//     for (int i = 0; i < set->size; i++)
//     {
//         int map0idx    = ((int *)set->cell_index_dat->data)[i * set->cell_index_dat->dim];

//         const int map1idx = arg1.map_data[map0idx * arg1.dim + 0];
//         const int map2idx = arg1.map_data[map0idx * arg1.dim + 1];
//         const int map3idx = arg1.map_data[map0idx * arg1.dim + 2];
//         const int map4idx = arg1.map_data[map0idx * arg1.dim + 3];

//         const int map5idx = arg5.map_data[map0idx * arg5.dim + 0];
//         const int map6idx = arg5.map_data[map0idx * arg5.dim + 1];
//         const int map7idx = arg5.map_data[map0idx * arg5.dim + 2];
//         const int map8idx = arg5.map_data[map0idx * arg5.dim + 3];

//         weight_particle_to_mesh_nodes__kernel(
//             &((double*)arg0.data)[i * arg0.dim],    // part_lc
//             &((double*)arg1.data)[map1idx],         // node_charge_den0
//             &((double*)arg1.data)[map2idx],         // node_charge_den1
//             &((double*)arg1.data)[map3idx],         // node_charge_den2
//             &((double*)arg1.data)[map4idx],         // node_charge_den3
//             &((double*)arg5.data)[map5idx],         // node_volume0
//             &((double*)arg5.data)[map6idx],         // node_volume1
//             &((double*)arg5.data)[map7idx],         // node_volume2
//             &((double*)arg5.data)[map8idx]          // node_volume3
//         );
//     }
// }

//*************************************************************************************************