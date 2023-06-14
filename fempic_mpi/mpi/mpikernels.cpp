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

//****************************************
double CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0;

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_spwt"))              CONST_spwt = *((double*)data);
    else if (!strcmp(name,"CONST_ion_velocity")) CONST_ion_velocity = *((double*)data);
    else if (!strcmp(name,"CONST_dt"))           CONST_dt = *((double*)data);
    else if (!strcmp(name,"CONST_plasma_den"))   CONST_plasma_den = *((double*)data);
    else if (!strcmp(name,"CONST_mass"))         CONST_mass = *((double*)data);
    else if (!strcmp(name,"CONST_charge"))       CONST_charge = *((double*)data);
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

//*************************************************************************************************
void opp_loop_inject__InjectIons(
    opp_set set,      // particles_set
    opp_arg arg0,     // part_position,
    opp_arg arg1,     // part_velocity,
    opp_arg arg2,     // part_cell_connectivity,
    opp_arg arg3,     // iface to cell map
    opp_arg arg4,     // cell_ef,
    opp_arg arg5,     // iface_u,
    opp_arg arg6,     // iface_v,
    opp_arg arg7,     // iface_normal,
    opp_arg arg8,     // iface_node_pos
    opp_arg arg9      // dummy_part_random
)
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_inject__InjectIons set_size %d diff %d", 
        set->size, set->diff);
    
    opp_profiler->start("InjectIons");

    const int inj_start = (set->size - set->diff);

    const int nargs = 10;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);
    args[6] = std::move(arg6);
    args[7] = std::move(arg7);
    args[8] = std::move(arg8);
    args[9] = std::move(arg9);

    int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args); // 

    int map0idx = -1, map1idx = 0, map0idxBackup = 0;
    int counter = 0;

    if (set_size > 0) 
    {
        for (int i = 0; i < set->diff; i++)
        {    
            map0idx = ((int *)set->mesh_relation_dat->data)[inj_start + i]; // iface index
            map1idx = args[4].map_data[map0idx]; // cell index

            // this is used to get the random numbers from the begining if the iface/cell change
            if (map0idx != map0idxBackup) 
            {
                map0idxBackup = map0idx;
                counter = 0;
            }

// {
//  opp_printf("general", "i=%d mesh mapping=%d part=%d", i, map0idx, (inj_start+i));
//  opp_printf("dummy_part_random", "%+2.20lE %+2.20lE", 
//       ((double*) args[9].data)[i * args[9].dim], ((double*) args[9].data)[i * args[9].dim+1]);
//  opp_printf("iface_node_pos", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[8].data)[map0idx * args[8].dim], 
//       ((double*) args[8].data)[map0idx * args[8].dim+1], ((double*) args[8].data)[map0idx * args[8].dim+2]);
//  opp_printf("iface_u", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[5].data)[map0idx * args[5].dim], 
//         ((double*) args[5].data)[map0idx * args[5].dim+1], ((double*) args[5].data)[map0idx * args[5].dim+2]);
//  opp_printf("iface_v", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[6].data)[map0idx * args[6].dim], 
//         ((double*) args[6].data)[map0idx * args[6].dim+1], ((double*) args[6].data)[map0idx * args[6].dim+2]);
//  opp_printf("iface_normal", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[7].data)[map0idx * args[7].dim], 
//         ((double*) args[7].data)[map0idx * args[7].dim+1], ((double*) args[7].data)[map0idx * args[7].dim+2]);
//  opp_printf("cell_ef", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[4].data)[map0idx * args[4].dim], 
//         ((double*) args[4].data)[map0idx * args[4].dim+1], ((double*) args[4].data)[map0idx * args[4].dim+2]);
// }
            inject_ions__kernel(
                &((double*) args[0].data)[(inj_start + i) * args[0].dim],     // part_position,
                &((double*) args[1].data)[(inj_start + i) * args[1].dim],     // part_velocity,
                &((int*)    args[2].data)[(inj_start + i) * args[2].dim],     // part_cell_connectivity,
                &((int*)    args[3].data)[map0idx * args[3].dim],             // iface to cell map
                &((double*) args[4].data)[map1idx * args[4].dim],             // cell_ef,
                &((double*) args[5].data)[map0idx * args[5].dim],             // iface_u,
                &((double*) args[6].data)[map0idx * args[6].dim],             // iface_v,
                &((double*) args[7].data)[map0idx * args[7].dim],             // iface_normal,
                &((double*) args[8].data)[map0idx * args[8].dim],             // iface_node_pos
                &((double*) args[9].data)[(counter++) * args[9].dim]          // dummy_part_random
            );
        }
    }

    opp_mpi_set_dirtybit(nargs, args);

    opp_profiler->end("InjectIons");
}

// #define DEBUG_INTERNAL

//*************************************************************************************************
void opp_loop_all_part_move__MoveToCells(
    opp_set set,      // particles_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // part_pos,
    opp_arg arg2,     // part_vel,
    opp_arg arg3,     // part_lc,
    opp_arg arg4,     // current_cell_index,
    opp_arg arg5,     // current_cell_volume,
    opp_arg arg6,     // current_cell_det,
    opp_arg arg7,     // cell_connectivity,
    opp_arg arg8,     // node_charge_den0,
    opp_arg arg9,     // node_charge_den1,
    opp_arg arg10,    // node_charge_den2,
    opp_arg arg11     // node_charge_den3,
)
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all_part_move__MoveToCells set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("MoveToCells");

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    double kernel_time = 0.0;
    int total_particles = 0;
    int comm_iteration = 0;
    int max_internal_hops = 0, internal_hops = 0;
    std::vector<int> particle_loops_per_comm_iter(10, 0);
    auto total_start = std::chrono::system_clock::now();
    auto start = std::chrono::system_clock::now();
#endif // -------------------------------------------------------------------------------------------

    int nargs = 12;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);
    args[5]  = std::move(arg5);
    args[6]  = std::move(arg6);
    args[7]  = std::move(arg7);
    args[8]  = std::move(arg8);
    args[9]  = std::move(arg9);
    args[10] = std::move(arg10);
    args[11] = std::move(arg11);

    int *map0idx = nullptr;

    int set_size = opp_mpi_halo_exchanges(set, nargs, args);

    // if there is access to a dat with OPP_Map_from_Mesh_Rel and a mapping, 
    // then we should reduce the contributions to the element containing rank
    // Here we should make the values of that dat to zero prior loop, 
    // execute the loop and communicate the outcome to the residing rank, like in a halo exchange, 
    // but when received, that rank should do the reduction
    opp_init_double_indirect_reductions(nargs, args);

    // unable to overlap much of computation and communication
    opp_mpi_halo_wait_all(nargs, args); 

    do // iterate until all mpi ranks say, I am done
    {
        opp_init_particle_move(set, nargs, args);
        
        if (FP_DEBUG) 
            opp_printf("FEMPIC", "opp_loop_particle_all__MoveToCells Starting iteration %d, start[%d] end[%d]", 
                OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        start = std::chrono::system_clock::now();
#endif // -------------------------------------------------------------------------------------------

        for (int i = OPP_iter_start; i < OPP_iter_end; i++)
        {
            opp_move_var m = opp_get_move_var();

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
            internal_hops = 0;
#endif // -------------------------------------------------------------------------------------------

            do
            { 
                map0idx = &(OPP_mesh_relation_data[i]);

                const int map1idx = args[8].map_data[*map0idx * args[8].map->dim + 0];
                const int map2idx = args[8].map_data[*map0idx * args[8].map->dim + 1];
                const int map3idx = args[8].map_data[*map0idx * args[8].map->dim + 2];
                const int map4idx = args[8].map_data[*map0idx * args[8].map->dim + 3];

                move_all_particles_to_cell__kernel(
                    (m),
                    &((double*) args[0].data)[*map0idx * args[0].dim],  // cell_ef,
                    &((double*) args[1].data)[i * args[1].dim],         // part_pos,
                    &((double*) args[2].data)[i * args[2].dim],         // part_vel,
                    &((double*) args[3].data)[i * args[3].dim],         // part_lc,
                    &((int *)   args[4].data)[i * args[4].dim],         // current_cell_index,
                    &((double*) args[5].data)[*map0idx * args[5].dim],  // current_cell_volume,
                    &((double*) args[6].data)[*map0idx * args[6].dim],  // current_cell_det,
                    &((int*)    args[7].data)[*map0idx * args[7].dim],  // cell_connectivity,
                    &((double*) args[8].data)[map1idx],                 // node_charge_den0,
                    &((double*) args[8].data)[map2idx],                 // node_charge_den1,
                    &((double*) args[8].data)[map3idx],                 // node_charge_den2,
                    &((double*) args[8].data)[map4idx]                  // node_charge_den3,
                );

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
                internal_hops++; // can remove
#endif // -------------------------------------------------------------------------------------------

            // opp_part_check_status checks whether map0idx is in halo list, 
            // if yes, pack the particle into MPI buffer and set status to NEED_REMOVE
            } while (opp_part_check_status(m, *map0idx, set, i, set->particle_remove_count));

#ifdef DEBUG_INTERNAL
            if (max_internal_hops < internal_hops) 
                max_internal_hops = internal_hops; // can remove
#endif // -------------------------------------------------------------------------------------------
        }

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
        kernel_time += (double)diff.count();
        total_particles += (OPP_iter_end - OPP_iter_start);
        particle_loops_per_comm_iter[OPP_comm_iteration] = (OPP_iter_end - OPP_iter_start);
        comm_iteration++;
#endif // -------------------------------------------------------------------------------------------

    } while (opp_finalize_particle_move(set)); // iterate until all mpi ranks say, I am done

    opp_exchange_double_indirect_reductions(nargs, args);

    // if auto_sort is set, then could sort here

    // TODO : can this be added to opp_mpi_halo_exchanges, to complete before the next usage?
    opp_complete_double_indirect_reductions(nargs, args);

    opp_mpi_set_dirtybit(nargs, args);

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    std::chrono::duration<double> total_diff   = std::chrono::system_clock::now() - total_start;
    opp_printf("MoveToCells", "TotalTime: %2.15lE KernelTime: %2.15lE | total_particles: %d | \
        particle_loops_per_comm_iter [%d %d %d %d] | comm_iteration: %d max_internal_hops: %d", 
        (double)total_diff.count(), kernel_time, total_particles, 
        particle_loops_per_comm_iter[0], particle_loops_per_comm_iter[1], particle_loops_per_comm_iter[2], 
        particle_loops_per_comm_iter[3], comm_iteration, max_internal_hops);
#endif // -------------------------------------------------------------------------------------------

    opp_profiler->end("MoveToCells");
}

//*************************************************************************************************
void opp_loop_all__ComputeNodeChargeDensity(
    opp_set set,     // nodes_set
    opp_arg arg0,    // node_charge_density
    opp_arg arg1     // node_volume
)
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__ComputeNodeChargeDensity set_size %d", set->size);

    opp_profiler->start("ComputeNodeChargeDensity");

    int nargs = 2;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);

    int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    if (set_size > 0) 
    {
        for (int i=0; i < set_size; i++) 
        {
            if (i == set->core_size) 
            {
                opp_mpi_halo_wait_all(nargs, args);
            }

            compute_node_charge_density__kernel(
                &((double*)arg0.data)[i * arg0.dim],
                &((double*)arg1.data)[i * arg1.dim]
            );
        }
    }  

    opp_mpi_set_dirtybit(nargs, args);

    opp_profiler->end("ComputeNodeChargeDensity");
}

//*************************************************************************************************
void opp_loop_all__ComputeElectricField(
    opp_set set,      // cells_set
    opp_arg arg0,     // cell_electric_field,
    opp_arg arg1,     // cell_shape_deriv,
    opp_arg arg2,     // node_potential0,
    opp_arg arg3,     // node_potential1,
    opp_arg arg4,     // node_potential2,
    opp_arg arg5      // node_potential3,
)
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__ComputeElectricField set_size %d", set->size);

    opp_profiler->start("ComputeElectricField");

    int nargs = 6;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);

    int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    if (set_size > 0) 
    {
        for (int i = 0; i < set_size; i++)
        {
            if (i == set->core_size) 
            {
                opp_mpi_halo_wait_all(nargs, args);
            }

            const int map1idx = arg2.map_data[i * arg2.map->dim + 0];
            const int map2idx = arg2.map_data[i * arg2.map->dim + 1];
            const int map3idx = arg2.map_data[i * arg2.map->dim + 2];
            const int map4idx = arg2.map_data[i * arg2.map->dim + 3];

            compute_electric_field__kernel(
                &((double*)arg0.data)[i * arg0.dim],    // cell_electric_field
                &((double*)arg1.data)[i * arg1.dim],    // cell_shape_deriv
                &((double*)arg2.data)[map1idx],         // node_potential0
                &((double*)arg2.data)[map2idx],         // node_potential1
                &((double*)arg2.data)[map3idx],         // node_potential2
                &((double*)arg2.data)[map4idx]          // node_potential3
            );
        } 
    }  

    opp_mpi_set_dirtybit(nargs, args);

    opp_profiler->end("ComputeElectricField");
}

//*************************************************************************************************