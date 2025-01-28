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

#include <opp_sycl.h>

OPP_INT* OPP_cells_set_size_d = nullptr;
OPP_INT OPP_cells_set_size = -1;

OPP_INT* OPP_comm_iteration_d = nullptr;
OPP_INT OPP_comm_iteration_h = -1;

OPP_INT* cellMapper_pos_stride_d = nullptr;
OPP_INT cellMapper_pos_stride = -1;

OPP_INT* OPP_rank_d;
OPP_INT* OPP_comm_size_d;

size_t* opp_maxSavedDHGrid_d = nullptr;
size_t* opp_minSavedDHGrid_d = nullptr;

// bool initialized = false;
// void opp_register_constants() {
//     if (!initialized) {
//         opp_register_const<OPP_INT>(cellMapper_pos_stride_d, 1);
//         opp_register_const<OPP_INT>(OPP_rank_d, 1);
//         opp_register_const<OPP_REAL>(OPP_comm_size_d, 1);
//         opp_register_const<size_t>(opp_maxSavedDHGrid_d, 3);
//         opp_register_const<size_t>(opp_minSavedDHGrid_d, 3);

//         initialized = true;
//     }
// }

//*******************************************************************************
// Returns true only if another hop is required by the current rank
inline bool opp_part_check_status_device(char& move_flag, bool& iter_one_flag, 
        int* cell_id, const int particle_index, int* remove_count, int *remove_part_indices, 
        int *move_part_indices, int *move_cell_indices, int *move_count, const OPP_INT* cell_set_size_d)
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE) {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE) {
        cell_id[0] = MAX_CELL_INDEX;
        const int removeArrayIndex = opp_atomic_fetch_add(remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (cell_id[0] >= cell_set_size_d[0])
    {
        // cell_id cell is not owned by the current mpi rank, need to communicate
        const int moveArrayIndex = opp_atomic_fetch_add(move_count, 1);
        move_part_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *cell_id;

        // To be removed from the current rank, packing will be done prior exchange & removal
        move_flag = OPP_NEED_REMOVE; // not required
        const int removeArrayIndex = opp_atomic_fetch_add(remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // cell_id is an own cell and move_flag == OPP_NEED_MOVE
    return true;
}

//*******************************************************************************
template <opp_access reduction, int intel, class T, class out_acc, class local_acc>
void opp_reduction(out_acc dat_g, int offset, T dat_l, local_acc temp, sycl::nd_item<1> &item) {
    T dat_t;

    /* important to finish all previous activity */
    item.barrier(sycl::access::fence_space::local_space); 

    size_t tid = item.get_local_id(0);
    temp[tid] = dat_l;

    for (size_t d = item.get_local_range(0) / 2; d > 0; d >>= 1) {
        item.barrier(sycl::access::fence_space::local_space);
        if (tid < d) {
        dat_t = temp[tid + d];

        switch (reduction) {
        case OPP_INC:
            dat_l = dat_l + dat_t;
            break;
        case OPP_MIN:
            if (dat_t < dat_l)
                dat_l = dat_t;
            break;
        case OPP_MAX:
            if (dat_t > dat_l)
                dat_l = dat_t;
            break;
        }
        temp[tid] = dat_l;
        }
    }

    if (tid == 0) {
        switch (reduction) {
        case OPP_INC:
            dat_g[offset] = dat_g[offset] + dat_l;
            break;
        case OPP_MIN:
            if (dat_l < dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        case OPP_MAX:
            if (dat_l > dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        }
    }
}
