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


#include "oppic_seq.h"
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
    int map0idx = -1, map1idx = 0, map0idxBackup = 0;
    int counter = 0;

    for (int i = 0; i < set->diff; i++)
    {    
        map0idx = ((int *)set->mesh_relation_dat->data)[inj_start + i]; // iface index
        map1idx = arg4.map_data[map0idx]; // cell index

        // this is used to get the random numbers from the begining if the iface/cell change
        if (map0idx != map0idxBackup) 
        {
            map0idxBackup = map0idx;
            counter = 0;
        }

        inject_ions__kernel(
            &((double *)arg0.data)[(inj_start + i) * arg0.dim],    // part_position,
            &((double *)arg1.data)[(inj_start + i) * arg1.dim],    // part_velocity,
            &((int *)arg2.data)[(inj_start + i) * arg2.dim],       // part_cell_connectivity,
            &((int *)arg3.data)[map0idx * arg3.dim],               // iface to cell map
            &((double*)arg4.data)[map1idx * arg4.dim],             // cell_ef,
            &((double*)arg5.data)[map0idx * arg5.dim],             // iface_u,
            &((double*)arg6.data)[map0idx * arg6.dim],             // iface_v,
            &((double*)arg7.data)[map0idx * arg7.dim],             // iface_normal,
            &((double*)arg8.data)[map0idx * arg8.dim],             // iface_node_pos
            &((double*)arg9.data)[(counter++) * arg9.dim]          // dummy_part_random
        );
    }

    opp_profiler->end("InjectIons");
}

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

    int *map0idx = nullptr;

    opp_init_particle_move(set, 0, nullptr);

    for (int i = 0; i < set->size; i++)
    {        
        opp_move_var m = opp_get_move_var();

        do
        { 
            m.OPP_inside_cell = true;

            map0idx = &(OPP_mesh_relation_data[i]);

            const int map1idx = arg8.map_data[*map0idx * arg8.map->dim + 0];
            const int map2idx = arg8.map_data[*map0idx * arg8.map->dim + 1];
            const int map3idx = arg8.map_data[*map0idx * arg8.map->dim + 2];
            const int map4idx = arg8.map_data[*map0idx * arg8.map->dim + 3];

            move_all_particles_to_cell__kernel(
                (m),
                &((double*) arg0.data)[*map0idx * arg0.dim],  // cell_ef,
                &((double*) arg1.data)[i * arg1.dim],         // part_pos,
                &((double*) arg2.data)[i * arg2.dim],         // part_vel,
                &((double*) arg3.data)[i * arg3.dim],         // part_lc,
                &((int *)   arg4.data)[i * arg4.dim],         // current_cell_index,
                &((double*) arg5.data)[*map0idx * arg5.dim],  // current_cell_volume,
                &((double*) arg6.data)[*map0idx * arg6.dim],  // current_cell_det,
                &((int*)    arg7.data)[*map0idx * arg7.dim],  // cell_connectivity,
                &((double*) arg8.data)[map1idx],                 // node_charge_den0,
                &((double*) arg8.data)[map2idx],                 // node_charge_den1,
                &((double*) arg8.data)[map3idx],                 // node_charge_den2,
                &((double*) arg8.data)[map4idx]                  // node_charge_den3,
            );         
            
            m.OPP_iteration_one = false;

        } while (opp_part_check_status(m, *map0idx, set, i, set->particle_remove_count));
    }

    opp_finalize_particle_move(set);

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

    for (int i=0; i<set->size; i++) 
    {
        compute_node_charge_density__kernel(
            &((double*)arg0.data)[i * arg0.dim],
            &((double*)arg1.data)[i * arg1.dim]
        );
    }

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

    for (int i = 0; i < set->size; i++)
    {
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

    opp_profiler->end("ComputeElectricField");   
}
