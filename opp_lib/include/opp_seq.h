
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

#pragma once

#include <opp_lib.h>


template <typename... T, typename... OPARG>
void opp_par_loop(void (*kernel)(T *...), char const *name, opp_set set, opp_iterate_type iter_type,
                 OPARG... arguments) {
    printf("opp_par_loop %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}

template <typename... T, typename... OPARG>
void opp_par_loop_particle(void (*kernel)(T *...), char const *name, opp_set set, opp_iterate_type iter_type,
                 OPARG... arguments) {
    printf("opp_par_loopp_particle %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}


inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
    (void)args;
    (void)data;
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
    (void)args;
    (void)data;
}

enum oppx_move_status : char
{
    OPPX_MOVE_DONE = 0,
    OPPX_NEED_MOVE,
    OPPX_NEED_REMOVE,
};

extern char opp_move_status_flag;
extern bool opp_move_hop_iter_one_flag;

#define OPP_PARTICLE_MOVE_DONE { opp_move_status_flag = OPPX_MOVE_DONE; }
#define OPP_PARTICLE_NEED_MOVE { opp_move_status_flag = OPPX_NEED_MOVE; }
#define OPP_PARTICLE_NEED_REMOVE { opp_move_status_flag = OPPX_NEED_REMOVE; }
#define OPP_DO_ONCE (opp_move_hop_iter_one_flag)
#define OPP_MOVE_RESET_FLAGS { opp_move_status_flag = OPPX_MOVE_DONE; opp_move_hop_iter_one_flag = true; }

//*************************************************************************************************
inline bool opp_check_part_move_status(const OPP_INT map0idx, const OPP_INT particle_index, int& remove_count) 
{
    opp_move_hop_iter_one_flag = false;

    if (opp_move_status_flag == OPPX_MOVE_DONE)
    {
        return false;
    }
    else if (opp_move_status_flag == OPPX_NEED_REMOVE)
    {
        remove_count += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }

    // map0idx is an own cell and opp_move_status_flag == OPPX_NEED_MOVE
    return true;
}


//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove(opp_set set, const opp_point& point, const int partIndex, int& cellIdx) {

    // Since SEQ use global indices, we can simply use findClosestGlobalCellIndex
    const size_t structCellIdx = cellMapper->findStructuredCellIndex(point);
    
    if (structCellIdx == MAX_CELL_INDEX) { // This happens when point is out of the unstructured mesh
        if (OPP_DBG)
            opp_printf("move", 
            "Remove %d [Struct cell index invalid - strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                partIndex, structCellIdx, point.x, point.y, point.z);

        cellIdx = MAX_CELL_INDEX;
        return true;
    }

    cellIdx = cellMapper->findClosestCellIndex(structCellIdx);           
    if (cellIdx == MAX_CELL_INDEX) { // Particle is outside the mesh, need to remove
        return true;
    }

    return false;
}