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


#include "opp_mpi.h"
#include "../fempic.h"

//****************************************
double CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0;
void oppic_decl_const_impl(int dim, int size, char* data, const char* name)
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
void oppic_inject__Increase_particle_count
(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
)
{ TRACE_ME;

    if (FP_DEBUG) opp_printf("FEMPIC", "oppic_inject__Increase_particle_count set_size %d diff %d", set->size, set->diff);

    for (int i = 0; i < set->size; i++)
    {   
        calculate_injection_distribution(
            ((int *)arg0.data),
            &((double *)arg1.data)[i],
            &((int *)arg2.data)[i],
            ((double *)arg3.data) 
        );
    }

    oppic_increase_particle_count(particles_set, *((int *)arg0.data));

    int* part_mesh_connectivity = (int *)particles_set->mesh_relation_dat->data;
    int* distribution           = (int *)arg2.data;

    int start = (particles_set->size - particles_set->diff);
    int j = 0;

    for (int i = 0; i < particles_set->diff; i++)
    {
        if (i >= distribution[j]) j++; // check whether it is j or j-1    
        part_mesh_connectivity[start + i] = j;
    }   
}

//*************************************************************************************************
void oppic_par_loop_inject__InjectIons(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_velocity,
    oppic_arg arg2,     // part_cell_connectivity,
    oppic_arg arg3,     // iface to cell map
    oppic_arg arg4,     // cell_ef,
    oppic_arg arg5,     // iface_u,
    oppic_arg arg6,     // iface_v,
    oppic_arg arg7,     // iface_normal,
    oppic_arg arg8,     // iface_node_pos
    oppic_arg arg9      // dummy_part_random
)
{ TRACE_ME;

    if (FP_DEBUG) opp_printf("FEMPIC", "oppic_par_loop_inject__InjectIons set_size %d diff %d", set->size, set->diff);

    int inj_start = (set->size - set->diff);

    int nargs = 10;
    oppic_arg args[nargs];

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

    int set_size = oppic_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args); // 

    int map0idx = -1, map0idxBackup = 0;
    int counter = 0;

    if (set_size > 0) 
    {
        for (int i = 0; i < set->diff; i++)
        {    
            map0idx    = ((int *)set->mesh_relation_dat->data)[inj_start + i]; // iface index
            int map1idx = args[4].map_data[map0idx]; // cell index

            if (map0idx != map0idxBackup)
            {
                // opp_printf("InjectIons", "making counter to zero old %d new %d counter %d", map0idxBackup, map0idx, counter);
                map0idxBackup = map0idx;
                counter = 0;
            }
// if (i == 0)
// {
//     opp_printf("general", "i=%d mesh mapping=%d part=%d", i, map0idx, (inj_start + i));
    
//     opp_printf("dummy_part_random", "%+2.20lE %+2.20lE", ((double*) args[9].data)[i * args[9].dim], ((double*) args[9].data)[i * args[9].dim + 1]);
//     opp_printf("iface_node_pos", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[8].data)[map0idx * args[8].dim], ((double*) args[8].data)[map0idx * args[8].dim + 1], ((double*) args[8].data)[map0idx * args[8].dim + 2]);
//     opp_printf("iface_u", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[5].data)[map0idx * args[5].dim], ((double*) args[5].data)[map0idx * args[5].dim + 1], ((double*) args[5].data)[map0idx * args[5].dim + 2]);
//     opp_printf("iface_v", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[6].data)[map0idx * args[6].dim], ((double*) args[6].data)[map0idx * args[6].dim + 1], ((double*) args[6].data)[map0idx * args[6].dim + 2]);
//     opp_printf("iface_normal", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[7].data)[map0idx * args[7].dim], ((double*) args[7].data)[map0idx * args[7].dim + 1], ((double*) args[7].data)[map0idx * args[7].dim + 2]);
//     opp_printf("cell_ef", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[4].data)[map0idx * args[4].dim], ((double*) args[4].data)[map0idx * args[4].dim + 1], ((double*) args[4].data)[map0idx * args[4].dim + 2]);
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
                &((double*) args[9].data)[(counter++) * args[9].dim]                    // dummy_part_random
            );

// if (i == 0)
// {
//     opp_printf("part_position", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[0].data)[(inj_start + i) * args[0].dim], ((double*) args[0].data)[(inj_start + i) * args[0].dim+1], ((double*) args[0].data)[(inj_start + i) * args[0].dim+2]);
//     opp_printf("part_velocity", "%+2.20lE %+2.20lE %+2.20lE", ((double*) args[1].data)[(inj_start + i) * args[1].dim], ((double*) args[1].data)[(inj_start + i) * args[1].dim+1], ((double*) args[1].data)[(inj_start + i) * args[1].dim+2]);
//     opp_printf("part_cell_connectivity", "%d", ((int*)    args[2].data)[(inj_start + i) * args[2].dim]);
// }

        }
    }

    oppic_mpi_set_dirtybit(nargs, args);
}

//*************************************************************************************************
void oppic_par_loop_particle_all__MoveToCells(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // cell_ef,
    oppic_arg arg1,     // part_pos,
    oppic_arg arg2,     // part_vel,
    oppic_arg arg3,     // part_lc,
    oppic_arg arg4,     // current_cell_index,
    oppic_arg arg5,     // current_cell_volume,
    oppic_arg arg6,     // current_cell_det,
    oppic_arg arg7,     // cell_connectivity,
    oppic_arg arg8,     // node_charge_den0,
    oppic_arg arg9,     // node_charge_den1,
    oppic_arg arg10,    // node_charge_den2,
    oppic_arg arg11     // node_charge_den3,
)
{ TRACE_ME;

    if (FP_DEBUG) opp_printf("FEMPIC", "oppic_par_loop_particle_all__MoveToCells set_size %d diff %d", set->size, set->diff);

    int nargs = 12;
    oppic_arg args[nargs];

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

// if there is access to a dat with OPP_Map_from_Mesh_Rel and a mapping, 
// then we should reduce the contributions to the element containing rank
// Here we should make the values of that dat to zero prior loop, 
// execute the loop and communicate the outcome to the residing rank, like in a halo exchange, 
// but when received, that rank should do the reduction

    //if (false) // incase if we want to print current ranks node_charge_density including halos
    // {
    //     std::string f = std::string("MoveToCell_B_") + std::to_string(OPP_my_rank);
    //     oppic_print_dat_to_txtfile(args[8].dat  , f.c_str(), "node_charge_density.dat");
    // }

    opp_init_double_indirect_reductions(nargs, args);

    //if (false) // incase if we want to print current ranks node_charge_density including halos
    // {
    //     std::string f = std::string("MoveToCell_A_") + std::to_string(OPP_my_rank);
    //     oppic_print_dat_to_txtfile(args[8].dat  , f.c_str(), "node_charge_density.dat");
    // }

    int set_size = oppic_mpi_halo_exchanges(set, nargs, args);

    // unable to overlap computation and communication, could overlap if particles are sorted according to cell index
    opp_mpi_halo_wait_all(nargs, args); 
int comm_iteration = 0;
    do // iterate until all mpi ranks say, I am done
    {
        oppic_init_particle_move(set, nargs, args);
        
        if (FP_DEBUG) opp_printf("FEMPIC", "oppic_par_loop_particle_all__MoveToCells Starting iteration %d, start[%d] end[%d]", 
            OPP_comm_iteration, OPP_iter_start, OPP_iter_end);
            
        for (int i = OPP_iter_start; i < OPP_iter_end; i++)
        {
            opp_move_var m = opp_get_move_var();

            do
            { 
                map0idx = &(OPP_mesh_relation_data[i]);

                const int map1idx = args[8].map_data[*map0idx * args[8].map->dim + 0];
                const int map2idx = args[8].map_data[*map0idx * args[8].map->dim + 1];
                const int map3idx = args[8].map_data[*map0idx * args[8].map->dim + 2];
                const int map4idx = args[8].map_data[*map0idx * args[8].map->dim + 3];

                move_all_particles_to_cell__kernel(
                    (m),
                    &((double*) args[0].data)[*map0idx * args[0].dim],  // const double *cell_ef,
                    &((double*) args[1].data)[i * args[1].dim],         // double *part_pos,
                    &((double*) args[2].data)[i * args[2].dim],         // double *part_vel,
                    &((double*) args[3].data)[i * args[3].dim],         // double *part_lc,
                    &((int *)   args[4].data)[i * args[4].dim],         // int* current_cell_index,
                    &((double*) args[5].data)[*map0idx * args[5].dim],  // const double *current_cell_volume,
                    &((double*) args[6].data)[*map0idx * args[6].dim],  // const double *current_cell_det,
                    &((int*)    args[7].data)[*map0idx * args[7].dim],  // const int *cell_connectivity,
                    &((double*) args[8].data)[map1idx],                 // double *node_charge_den0,
                    &((double*) args[8].data)[map2idx],                 // double *node_charge_den1,
                    &((double*) args[8].data)[map3idx],                 // double *node_charge_den2,
                    &((double*) args[8].data)[map4idx]                  // double *node_charge_den3,
                );

                // should check whether map0idx is in halo list, if yes, pack the particle into MPI buffer and set status to NEED_REMOVE
            } while (opp_part_check_status(m, *map0idx, set, i, set->particle_remove_count));
        }
comm_iteration++;
    } while (oppic_finalize_particle_move(set)); // iterate until all mpi ranks say, I am done

if (OPP_my_rank == OPP_MPI_ROOT) // OPP_comm_iteration > 1 && 
    opp_printf("FEMPIC", "Multiple communication particle hops %d", comm_iteration);

    opp_exchange_double_indirect_reductions(nargs, args);

    // if auto_sort is set, then sort here

    // TODO : can this be added to oppic_mpi_halo_exchanges, to complete before the next usage?
    opp_complete_double_indirect_reductions(nargs, args);

    // TODO : Dirty bit should not be set for double indirect reductions, if opp_complete_double_indirect_reductions is called here
    oppic_mpi_set_dirtybit(nargs, args);
}

//*************************************************************************************************
void oppic_par_loop_all__ComputeNodeChargeDensity(
    oppic_set set,     // nodes_set
    oppic_arg arg0,    // node_charge_density
    oppic_arg arg1     // node_volume
)
{ TRACE_ME;
    
    if (FP_DEBUG) opp_printf("FEMPIC", "oppic_par_loop_all__ComputeNodeChargeDensity set_size %d", set->size);

    int nargs = 2;
    oppic_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    int set_size = oppic_mpi_halo_exchanges(set, nargs, args);
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

    oppic_mpi_set_dirtybit(nargs, args);
}

//*************************************************************************************************
void oppic_par_loop_all__ComputeElectricField(
    oppic_set set,      // cells_set
    oppic_arg arg0,     // cell_electric_field,
    oppic_arg arg1,     // cell_shape_deriv,
    oppic_arg arg2,     // node_potential0,
    oppic_arg arg3,     // node_potential1,
    oppic_arg arg4,     // node_potential2,
    oppic_arg arg5      // node_potential3,
)
{ TRACE_ME;

    if (FP_DEBUG) opp_printf("FEMPIC", "oppic_par_loop_all__ComputeElectricField set_size %d", set->size);

    int nargs = 6;
    oppic_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    int set_size = oppic_mpi_halo_exchanges(set, nargs, args);
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

    oppic_mpi_set_dirtybit(nargs, args);
}
