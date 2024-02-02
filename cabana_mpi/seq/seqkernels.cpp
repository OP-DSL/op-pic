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
#include "../cabana_defs.h"

const OPP_REAL CONST_qdt_2mc = -0.079733736813068389892578125;
const OPP_REAL CONST_cdt_dx = 0.637869894504547119140625;
const OPP_REAL CONST_cdt_dy = 0.757114291191101074218750;
const OPP_REAL CONST_cdt_dz = 0.15946747362613677978515625;
const OPP_REAL CONST_qsp = -1.0;

const OPP_REAL CONST_dx = +2.500000000000000000000000000000E-01;
const OPP_REAL CONST_dy = +2.106253922764312613846726662814E-01;
const OPP_REAL CONST_dz = +1.000000000000000000000000000000E+00;
const OPP_REAL CONST_dt = +1.594674816314074206058393201602E-01;

const OPP_REAL CONST_px = +6.378699265256296824233572806406E-01;
const OPP_REAL CONST_py = +7.571142297131837617385485827981E-01;
const OPP_REAL CONST_pz = +0.000000000000000000000000000000E+00;
const OPP_REAL CONST_cj = +1.594674736261367797851562500000E-01;

OPP_INT CONST_nx = -1;
OPP_INT CONST_ny = -1;
OPP_INT CONST_nz = -1;
OPP_INT CONST_ng = -1;

//****************************************
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_nx"))       CONST_nx = *((int*)data);
    else if (!strcmp(name,"CONST_ny"))  CONST_ny = *((int*)data);
    else if (!strcmp(name,"CONST_nz"))  CONST_nz = *((int*)data);
    else if (!strcmp(name,"CONST_ng"))  CONST_ng = *((int*)data);
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

//*************************************************************************************************
void opp_loop_all__interpolate_mesh_fields(
    opp_set set,        // cells_set
    opp_arg arg0,       // cell0_e,        // OPP_READ
    opp_arg arg1,       // cell0_b,        // OPP_READ
    opp_arg arg2,       // cell_x_e,       // OPP_READ
    opp_arg arg3,       // cell_y_e,       // OPP_READ
    opp_arg arg4,       // cell_z_e,       // OPP_READ
    opp_arg arg5,       // cell_yz_e,      // OPP_READ 
    opp_arg arg6,       // cell_xz_e,      // OPP_READ
    opp_arg arg7,       // cell_xy_e,      // OPP_READ
    opp_arg arg8,       // cell_x_b,       // OPP_READ
    opp_arg arg9,       // cell_y_b,       // OPP_READ
    opp_arg arg10,      // cell_z_b        // OPP_READ
    opp_arg arg11       // cell0_interp    // OPP_WRITE
)
{

    if (FP_DEBUG) opp_printf("CABANA", "opp_loop_all__interpolate_mesh_fields set_size %d", set->size);

    opp_profiler->start("Interpolate");

    for (int n = 0; n < set->size; n++)
    {
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::xu_y_z];
        const int map_3idx  = arg3.map_data[n * arg3.map->dim + CellMap::x_yu_z];
        const int map_4idx  = arg4.map_data[n * arg4.map->dim + CellMap::x_y_zu];
        const int map_5idx  = arg5.map_data[n * arg5.map->dim + CellMap::xu_yu_z];
        const int map_6idx  = arg6.map_data[n * arg6.map->dim + CellMap::x_yu_zu];
        const int map_7idx  = arg7.map_data[n * arg7.map->dim + CellMap::xu_y_zu];
        const int map_8idx  = arg8.map_data[n * arg8.map->dim + CellMap::xu_y_z];
        const int map_9idx  = arg9.map_data[n * arg9.map->dim + CellMap::x_yu_z];
        const int map_10idx = arg10.map_data[n * arg10.map->dim + CellMap::x_y_zu];

        interpolate_mesh_fields_kernel(
            &((double*) arg0.data)[n * arg0.dim],            // cell0_e,        // OPP_READ
            &((double*) arg1.data)[n * arg1.dim],            // cell0_b,        // OPP_READ
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_x_e,       // OPP_READ
            &((double*) arg3.data)[map_3idx * arg3.dim],     // cell_y_e,       // OPP_READ
            &((double*) arg4.data)[map_4idx * arg4.dim],     // cell_z_e,       // OPP_READ
            &((double*) arg5.data)[map_5idx * arg5.dim],     // cell_yz_e,      // OPP_READ
            &((double*) arg6.data)[map_6idx * arg6.dim],     // cell_xz_e,      // OPP_READ
            &((double*) arg7.data)[map_7idx * arg7.dim],     // cell_xy_e,      // OPP_READ
            &((double*) arg8.data)[map_8idx * arg8.dim],     // cell_x_b,       // OPP_READ
            &((double*) arg9.data)[map_9idx * arg9.dim],     // cell_y_b,       // OPP_READ
            &((double*) arg10.data)[map_10idx * arg10.dim],  // cell_z_b        // OPP_READ
            &((double*) arg11.data)[n * arg11.dim]           // cell0_interp    // OPP_WRITE
        );
    }

    opp_profiler->end("Interpolate");
}


//*************************************************************************************************
void opp_particle_mover__Move(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_cid           // OPP_RW
    opp_arg arg1,       // part_vel           // OPP_RW
    opp_arg arg2,       // part_pos           // OPP_RW
    opp_arg arg3,       // part_streak_mid    // OPP_RW
    opp_arg arg4,       // part_weight        // OPP_READ
    opp_arg arg5,       // cell_inter         // OPP_READ
    opp_arg arg6        // cell_acc           // OPP_INC
)
{

    if (FP_DEBUG) 
        opp_printf("CABANA", "opp_particle_mover__Move set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Move");

    opp_init_particle_move(set, 0, nullptr);

    int* cellIdx = nullptr;

    for (int n = 0; n < set->size; n++)
    {
        opp_move_var m; // = opp_get_move_var();

        do
        {
            cellIdx = &(OPP_mesh_relation_data[n]);

            push_particles_kernel(m, 
                &((int*)    arg0.data)[n * arg0.dim],        // part_cid 
                &((double*) arg1.data)[n * arg1.dim],        // part_vel 
                &((double*) arg2.data)[n * arg2.dim],        // part_pos 
                &((double*) arg3.data)[n * arg3.dim],        // part_streak_mid 
                &((double*) arg4.data)[n * arg4.dim],        // part_weight 
                &((double*) arg5.data)[*cellIdx * arg5.dim], // cell_interp 
                &((double*) arg6.data)[*cellIdx * arg6.dim]  // cell_acc
            );

        } while (opp_part_check_status(m, *cellIdx, set, n, set->particle_remove_count)); 
    }

    opp_finalize_particle_move(set);

    opp_profiler->end("Move");
}

//*************************************************************************************************
void opp_loop_all__accumulate_current_to_cells(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_j         // OPP_WRITE
    opp_arg arg1,    // cell0_acc       // OPP_READ
    opp_arg arg2,    // cell_xd_acc     // OPP_READ
    opp_arg arg3,    // cell_yd_acc     // OPP_READ
    opp_arg arg4,    // cell_zd_acc     // OPP_READ
    opp_arg arg5,    // cell_xyd_acc    // OPP_READ
    opp_arg arg6,    // cell_yzd_acc    // OPP_READ
    opp_arg arg7     // cell_xzd_acc    // OPP_READ
)
{

    if (FP_DEBUG) opp_printf("CABANA", "opp_loop_all__accumulate_current_to_cells set_size %d", set->size);

    opp_profiler->start("Acc_Current");

    for (int n = 0; n < set->size; n++)
    {
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::xd_y_z];
        const int map_3idx  = arg3.map_data[n * arg3.map->dim + CellMap::x_yd_z];
        const int map_4idx  = arg4.map_data[n * arg4.map->dim + CellMap::x_y_zd];
        const int map_5idx  = arg5.map_data[n * arg5.map->dim + CellMap::xd_yd_z];
        const int map_6idx  = arg6.map_data[n * arg6.map->dim + CellMap::x_yd_zd];
        const int map_7idx  = arg7.map_data[n * arg7.map->dim + CellMap::xd_y_zd];

        accumulate_current_to_cells_kernel(   
            &((double*) arg0.data)[n * arg0.dim],            // cell0_j     
            &((double*) arg1.data)[n * arg1.dim],            // cell0_acc   
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_xd_acc 
            &((double*) arg3.data)[map_3idx * arg3.dim],     // cell_yd_acc 
            &((double*) arg4.data)[map_4idx * arg4.dim],     // cell_zd_acc
            &((double*) arg5.data)[map_5idx * arg5.dim],     // cell_xyd_acc 
            &((double*) arg6.data)[map_6idx * arg6.dim],     // cell_yzd_acc 
            &((double*) arg7.data)[map_7idx * arg7.dim]      // cell_xzd_acc
        );
    }

    opp_profiler->end("Acc_Current");
}

//*************************************************************************************************
void opp_loop_all__half_advance_b(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_e        // OPP_WRITE
    opp_arg arg1,    // cell_y_e        // OPP_READ
    opp_arg arg2,    // cell_z_e        // OPP_READ
    opp_arg arg3,    // cell0_e         // OPP_READ
    opp_arg arg4     // cell0_b         // OPP_INC
)
{
    if (FP_DEBUG) opp_printf("CABANA", "opp_loop_all__half_advance_b set_size %d", set->size);

    opp_profiler->start("HalfAdv_B");

    for (int n = 0; n < set->size; n++)
    {
        const int map_0idx  = arg0.map_data[n * arg0.map->dim + CellMap::xu_y_z];
        const int map_1idx  = arg1.map_data[n * arg1.map->dim + CellMap::x_yu_z];
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::x_y_zu];

        half_advance_b_kernel (
            &((double*) arg0.data)[map_0idx * arg0.dim],     // cell_x_e, 
            &((double*) arg1.data)[map_1idx * arg1.dim],     // cell_y_e, 
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_z_e, 
            &((double*) arg3.data)[n * arg3.dim],            // cell0_e, 
            &((double*) arg4.data)[n * arg4.dim]             // cell0_b
        );
    }

    opp_profiler->end("HalfAdv_B");
}

//*************************************************************************************************
void opp_loop_all__advance_e(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_b        // OPP_READ
    opp_arg arg1,    // cell_y_b        // OPP_READ
    opp_arg arg2,    // cell_z_b        // OPP_READ
    opp_arg arg3,    // cell0_b         // OPP_READ
    opp_arg arg4,    // cell0_j         // OPP_READ
    opp_arg arg5     // cell0_e         // OPP_INC
)
{
    if (FP_DEBUG) opp_printf("CABANA", "opp_loop_all__advance_e set_size %d", set->size);

    opp_profiler->start("Adv_E");

    for (int n = 0; n < set->size; n++)
    {
        const int map_0idx  = arg0.map_data[n * arg0.map->dim + CellMap::xd_y_z];
        const int map_1idx  = arg1.map_data[n * arg1.map->dim + CellMap::x_yd_z];
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::x_y_zd];

        advance_e_kernel (
            &((double*) arg0.data)[map_0idx * arg0.dim],     // cell_x_b  
            &((double*) arg1.data)[map_1idx * arg1.dim],     // cell_y_b  
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_z_b  
            &((double*) arg3.data)[n * arg3.dim],            // cell0_b   
            &((double*) arg4.data)[n * arg4.dim],            // cell0_j   
            &((double*) arg5.data)[n * arg5.dim]             // cell0_e   
        );
    }

    opp_profiler->end("Adv_E");
}

