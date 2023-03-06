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

int moveToCells_all_stride_OPP_HOST_0 = -1;
int moveToCells_all_stride_OPP_HOST_1 = -1;
int moveToCells_all_stride_OPP_HOST_2 = -1;
int moveToCells_all_stride_OPP_HOST_3 = -1;
int moveToCells_all_stride_OPP_HOST_6 = -1;
int moveToCells_all_stride_OPP_HOST_7 = -1;
int moveToCells_all_stride_OPP_HOST_8 = -1;

__constant__ int moveToCells_all_stride_OPP_CUDA_0;
__constant__ int moveToCells_all_stride_OPP_CUDA_1;
__constant__ int moveToCells_all_stride_OPP_CUDA_2;
__constant__ int moveToCells_all_stride_OPP_CUDA_3;
__constant__ int moveToCells_all_stride_OPP_CUDA_6;
__constant__ int moveToCells_all_stride_OPP_CUDA_7;
__constant__ int moveToCells_all_stride_OPP_CUDA_8;

//user function
//*************************************************************************************************
__device__ void move_all_particles_to_cell__kernel(
    move_var* m,
    const double *cell_ef,
    double *part_pos,
    double *part_vel,
    double *part_lc,
    int* current_cell_index,
    const double *current_cell_volume,
    const double *current_cell_det,
    const int *cell_connectivity,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
)
{
    if (m->OPP_iteration_one)
    {
        for (int i = 0; i < DIMENSIONS; i++)
        {
            // double v0 = part_vel[i * moveToCells_all_stride_OPP_CUDA_2];
            // double v1 = (CONST_charge_cuda / CONST_mass_cuda);
            // double v2 = (cell_ef[i] * CONST_dt_cuda);
            // double v3 = (v1 * v2);
            // part_vel[i * moveToCells_all_stride_OPP_CUDA_2] = (v0 + v3);

            part_vel[i * moveToCells_all_stride_OPP_CUDA_2] += (CONST_charge_cuda / CONST_mass_cuda * cell_ef[i * moveToCells_all_stride_OPP_CUDA_0] * (CONST_dt_cuda));
// printf("%d %+2.20lE - %+2.20lE %+2.20lE\n - %+2.20lE %+2.20lE %+2.20lE %+2.20lE %+2.20lE\n\n", i,
//     part_vel[i * injectIons_stride_OPP_CUDA_1],
//     iface_normal[i * injectIons_stride_OPP_CUDA_7], cell_ef[i * injectIons_stride_OPP_CUDA_4],
//     v1, 
//     v2,
//     v3, 
//     v4,
//     v5);      
        
        }
   
        for (int i = 0; i < DIMENSIONS; i++)
            part_pos[i * moveToCells_all_stride_OPP_CUDA_1] += part_vel[i * moveToCells_all_stride_OPP_CUDA_2] * (CONST_dt_cuda); // v = u + at
    }

    for (int i=0; i<NODES_PER_CELL; i++) /*loop over vertices*/
    {
        part_lc[i * moveToCells_all_stride_OPP_CUDA_3] = (1.0/6.0) * (
            current_cell_det[(i * DET_FIELDS + 0) * moveToCells_all_stride_OPP_CUDA_6] - 
            current_cell_det[(i * DET_FIELDS + 1) * moveToCells_all_stride_OPP_CUDA_6] * part_pos[0 * moveToCells_all_stride_OPP_CUDA_1] + 
            current_cell_det[(i * DET_FIELDS + 2) * moveToCells_all_stride_OPP_CUDA_6] * part_pos[1 * moveToCells_all_stride_OPP_CUDA_1] - 
            current_cell_det[(i * DET_FIELDS + 3) * moveToCells_all_stride_OPP_CUDA_6] * part_pos[2 * moveToCells_all_stride_OPP_CUDA_1]
                ) / (*current_cell_volume);
        
        if (part_lc[i * moveToCells_all_stride_OPP_CUDA_3]<0 || part_lc[i * moveToCells_all_stride_OPP_CUDA_3]>1.0) m->OPP_inside_cell = false;
    }    

    if (m->OPP_inside_cell)
    {
        m->OPP_move_status = OPP_MOVE_DONE;

        atomicAdd(node_charge_den0, (part_lc[0 * moveToCells_all_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den1, (part_lc[1 * moveToCells_all_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den2, (part_lc[2 * moveToCells_all_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den3, (part_lc[3 * moveToCells_all_stride_OPP_CUDA_3]));

        return;
    }

    // outside the last known cell, find most negative weight and use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0 * moveToCells_all_stride_OPP_CUDA_3];
    
    for (int i=1; i<NEIGHBOUR_CELLS; i++)
    {
        if (part_lc[i * moveToCells_all_stride_OPP_CUDA_3] < min_lc) 
        {
            min_lc = part_lc[i * moveToCells_all_stride_OPP_CUDA_3];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i * moveToCells_all_stride_OPP_CUDA_7] >= 0) // is there a neighbor in this direction?
    {
        (*current_cell_index) = cell_connectivity[min_i * moveToCells_all_stride_OPP_CUDA_7];
        m->OPP_move_status = OPP_NEED_MOVE;
    }
    else
    {
        m->OPP_move_status = OPP_NEED_REMOVE;
    }
}

// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_all_MoveToCells(
    int *__restrict d_cell_index,
    double *__restrict ind_arg0,
    double *__restrict dir_arg1,
    double *__restrict dir_arg2,
    double *__restrict dir_arg3,
    int *__restrict dir_arg4,
    const double *__restrict ind_arg5,
    const double *__restrict ind_arg6,
    const int *__restrict ind_arg7,
    double *__restrict ind_arg8,
    const int *__restrict opDat8Map,
    double *__restrict ind_arg9,
    double *__restrict ind_arg10,
    double *__restrict ind_arg11,
    int *__restrict particle_remove_count,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;
        
        // int& map0idx = dir_arg4[n];
        // int& map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? dir_arg2 and d_cell_index has same pointer values, but this get stuck!

        move_var m;

        do
        {
            m.OPP_inside_cell = true;
            int map0idx = dir_arg4[n];

            const int map1idx = opDat8Map[map0idx + moveToCells_all_stride_OPP_CUDA_8 * 0];
            const int map2idx = opDat8Map[map0idx + moveToCells_all_stride_OPP_CUDA_8 * 1];
            const int map3idx = opDat8Map[map0idx + moveToCells_all_stride_OPP_CUDA_8 * 2];
            const int map4idx = opDat8Map[map0idx + moveToCells_all_stride_OPP_CUDA_8 * 3];

            move_all_particles_to_cell__kernel(
                &(m),
                (ind_arg0 + map0idx),   // cell_ef,
                (dir_arg1 + n),         // part_pos,
                (dir_arg2 + n),         // part_vel,
                (dir_arg3 + n),         // part_lc,
                (dir_arg4 + n),         // current_cell_index,
                (ind_arg5 + map0idx),   // current_cell_volume,
                (ind_arg6 + map0idx),   // current_cell_det,
                (ind_arg7 + map0idx),   // cell_connectivity,
                (ind_arg8 + map1idx),   // node_charge_den0,
                (ind_arg9 + map2idx),   // node_charge_den1,
                (ind_arg10 + map3idx),  // node_charge_den2,
                (ind_arg11 + map4idx)   // node_charge_den3,
            );                
            
            m.OPP_iteration_one = false;

// if (m.OPP_move_status != (int)OPP_NEED_MOVE)
// {
//     printf("%d %d - %d - %d %d %d %d -\n\t\t %+2.25lE %+2.25lE %+2.25lE %+2.25lE\n", n, map0idx, (int)m.OPP_move_status,
//         map1idx, map2idx, map3idx, map4idx,
//         (dir_arg3 + n)[0 * moveToCells_all_stride_OPP_CUDA_3],
//         (dir_arg3 + n)[1 * moveToCells_all_stride_OPP_CUDA_3],
//         (dir_arg3 + n)[2 * moveToCells_all_stride_OPP_CUDA_3],
//         (dir_arg3 + n)[3 * moveToCells_all_stride_OPP_CUDA_3]);                
// }

        } while (m.OPP_move_status == (int)OPP_NEED_MOVE);

        if (m.OPP_move_status == OPP_NEED_REMOVE) 
        {
            atomicAdd(particle_remove_count, 1);
            dir_arg4[n] = MAX_CELL_INDEX;
        }
    }
}

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
    
    if (FP_DEBUG) printf("FEMPIC - oppic_par_loop_particle_all__MoveToCells set_size %d\n", set->size);

    int nargs = 12;
    oppic_arg args[nargs] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11 };
    
    oppic_init_particle_move(set);

    int set_size = oppic_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        moveToCells_all_stride_OPP_HOST_0 = arg0.dat->set->set_capacity;
        moveToCells_all_stride_OPP_HOST_1 = arg1.dat->set->set_capacity;
        moveToCells_all_stride_OPP_HOST_2 = arg2.dat->set->set_capacity;
        moveToCells_all_stride_OPP_HOST_3 = arg3.dat->set->set_capacity;
        moveToCells_all_stride_OPP_HOST_6 = arg6.dat->set->set_capacity; 
        moveToCells_all_stride_OPP_HOST_7 = arg7.size;
        moveToCells_all_stride_OPP_HOST_8 = arg8.map->from->size;

        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_0, &moveToCells_all_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_1, &moveToCells_all_stride_OPP_HOST_1, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_2, &moveToCells_all_stride_OPP_HOST_2, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_3, &moveToCells_all_stride_OPP_HOST_3, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_6, &moveToCells_all_stride_OPP_HOST_6, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_7, &moveToCells_all_stride_OPP_HOST_7, sizeof(int));
        cudaMemcpyToSymbol(moveToCells_all_stride_OPP_CUDA_8, &moveToCells_all_stride_OPP_HOST_8, sizeof(int));

        int start   = 0;
        int end     = set->size;

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_all_MoveToCells<<<nblocks, nthread>>>(
                (int *)     set->mesh_relation_dat->data_d,
                (double *)  arg0.data_d,
                (double *)  arg1.data_d,
                (double *)  arg2.data_d,
                (double *)  arg3.data_d,
                (int *)     arg4.data_d,
                (double *)  arg5.data_d,
                (double *)  arg6.data_d,
                (int *)     arg7.data_d,
                (double *)  arg8.data_d,
                (int *)     arg8.map_data_d,
                (double *)  arg9.data_d,
                (double *)  arg10.data_d,
                (double *)  arg11.data_d,
                (int *)     set->particle_remove_count_d,
                start, 
                end);
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());

    oppic_finalize_particle_move(set);

    oppic_mpi_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());
}

//*************************************************************************************************