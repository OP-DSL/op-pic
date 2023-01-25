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

int opDat0_MoveToCells_inject_stride_OPPIC_HOST = -1;
int opDat1_MoveToCells_inject_stride_OPPIC_HOST = -1;
int opDat3_MoveToCells_inject_stride_OPPIC_HOST = -1;
int opDat5_MoveToCells_inject_stride_OPPIC_HOST = -1;
int opDat6_MoveToCells_inject_stride_OPPIC_HOST = -1;
int opDat7_MoveToCells_inject_stride_OPPIC_HOST = -1;

__constant__ int opDat0_MoveToCells_inject_stride_OPPIC_CONSTANT;
__constant__ int opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT;
__constant__ int opDat3_MoveToCells_inject_stride_OPPIC_CONSTANT;
__constant__ int opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT;
__constant__ int opDat6_MoveToCells_inject_stride_OPPIC_CONSTANT;
__constant__ int opDat7_MoveToCells_inject_stride_OPPIC_CONSTANT;


//user function
//*************************************************************************************************
__device__ void move_injected_particles_to_cell__kernel(
    int* move_status,
    const double* part_pos,
    double* part_lc,
    int* current_cell_index,
    double* part_vel,
    const double *cell_volume,
    const double *cell_det,
    const double *cell_ef,
    const int *cell_connectivity,
    const bool* search,
    const double* dt)
{
    bool inside = true;

    for (int i=0; i < NODES_PER_CELL; i++) /*loop over vertices*/
    {
        part_lc[i * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT] = (1.0/6.0) * (
            cell_det[(0 + i * NODES_PER_CELL) * opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT] - 
            cell_det[(1 + i * NODES_PER_CELL) * opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT] * part_pos[0 * opDat0_MoveToCells_inject_stride_OPPIC_CONSTANT] + 
            cell_det[(2 + i * NODES_PER_CELL) * opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT] * part_pos[1 * opDat0_MoveToCells_inject_stride_OPPIC_CONSTANT] - 
            cell_det[(3 + i * NODES_PER_CELL) * opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT] * part_pos[2 * opDat0_MoveToCells_inject_stride_OPPIC_CONSTANT]
                ) / (*cell_volume);
        
        if (part_lc[i * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT] < 0 || part_lc[i * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT] > 1.0) inside = false;
    }

    if (inside)
    {
        *move_status = MOVE_DONE;
        
        part_vel[0 * opDat3_MoveToCells_inject_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[0 * opDat6_MoveToCells_inject_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
    
        part_vel[1 * opDat3_MoveToCells_inject_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[1 * opDat6_MoveToCells_inject_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
    
        part_vel[2 * opDat3_MoveToCells_inject_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[2 * opDat6_MoveToCells_inject_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
        return;
    }

    if (*search) 
    {
        (*current_cell_index)++; // outside the last known cell, Increment the cell_index to search in the full mesh
        return;
    }

    // outside the last known cell, find most negative weight and use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0 * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT];
    
    for (int i=1; i<NEIGHBOUR_CELLS; i++)
    {
        if (part_lc[i * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT] < min_lc) 
        {
            min_lc = part_lc[i * opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i * opDat7_MoveToCells_inject_stride_OPPIC_CONSTANT] >= 0) // is there a neighbor in this direction?
    {
        (*current_cell_index) = cell_connectivity[min_i * opDat7_MoveToCells_inject_stride_OPPIC_CONSTANT];
        *move_status = NEED_MOVE;
    }
    else
    {
        (*current_cell_index) = MAX_CELL_INDEX;
        *move_status = NEED_REMOVE;
    }
}


// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_injected_MoveToCells(
    int *__restrict d_cell_index,
    double *__restrict dir_arg0,            // part_position,
    double *__restrict dir_arg1,            // part_weights,
    int *__restrict dir_arg2,               // part_cell_index,
    double *__restrict ind_arg3,            // part_velocity,
    const double *__restrict ind_arg4,      // cell_volume,
    const double *__restrict ind_arg5,      // cell_det,
    const double *__restrict ind_arg6,      // cell_electric_field,
    const int *__restrict ind_arg7,         // cell_connectivity_map,
    const bool *__restrict dir_arg8,        // particles_injected,
    const double *__restrict dir_arg9,      // dt,
    int *__restrict particle_statuses,      // mark particles as MOVE_DONE, NEED_REMOVE etc...
    int start,
    int end,
    int num_cells) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;
        
        int& map0idx = dir_arg2[n];
        // int& map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? dir_arg2 and d_cell_index has same pointer values, but this get stuck!

        int move_status = (int)NEED_MOVE;

        do
        {
            //user-supplied kernel call
            move_injected_particles_to_cell__kernel(
                &(move_status),
                (dir_arg0 + n),         // part_pos
                (dir_arg1 + n),         // part_weights
                (dir_arg2 + n),         // part_cell_index
                (ind_arg3 + n),         // part_velocity,
                (ind_arg4 + map0idx),   // cell_volume
                (ind_arg5 + map0idx),   // cell_det
                (ind_arg6 + map0idx),   // cell_electric_field
                (ind_arg7 + map0idx),   // cell_connectivity
                dir_arg8,               // full_mesh_search
                dir_arg9                // dt
            );                
            
        } while ((move_status == (int)NEED_MOVE) && (map0idx < num_cells));

        if (move_status == (int)NEED_REMOVE) /*outside the mesh*/
        {  
            printf("***** INJ : need to remove particle %d\n", n);            
            particle_statuses[n] = (int)NEED_REMOVE;
        }
        else if (move_status != MOVE_DONE)
        {
            printf("****************************** move_injected_particles_to_cell__kernel returned status[%d] for particle [%d]\n", move_status, n);
        }

        // __syncthreads();
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
    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_particle_inject__MoveToCells num_particles %d\n", set->size);

    int nargs = 10;
    oppic_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;
    args[7] = arg7;
    args[8] = arg8;
    args[9] = arg9;

    cutilSafeCall(cudaMalloc(&(arg8.data_d), arg8.size));
    cutilSafeCall(cudaMemcpy(arg8.data_d, arg8.data, arg8.size, cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMalloc(&(arg9.data_d), arg9.size));
    cutilSafeCall(cudaMemcpy(arg9.data_d, arg9.data, arg9.size, cudaMemcpyHostToDevice));

    oppic_init_particle_move(set);

    int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        opDat0_MoveToCells_inject_stride_OPPIC_HOST = arg0.dat->set->size;
        opDat1_MoveToCells_inject_stride_OPPIC_HOST = arg1.dat->set->size;
        opDat3_MoveToCells_inject_stride_OPPIC_HOST = arg3.dat->set->size;
        opDat5_MoveToCells_inject_stride_OPPIC_HOST = arg5.dat->set->size;
        opDat6_MoveToCells_inject_stride_OPPIC_HOST = arg6.dat->set->size;
        opDat7_MoveToCells_inject_stride_OPPIC_HOST = arg7.map->from->size;

        cudaMemcpyToSymbol(opDat0_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat0_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat1_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat1_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat3_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat3_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat5_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat5_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat6_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat6_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat7_MoveToCells_inject_stride_OPPIC_CONSTANT, &opDat7_MoveToCells_inject_stride_OPPIC_HOST, sizeof(int));

        int start   = (set->size - set->diff);
        int end     = set->size;

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_injected_MoveToCells<<<nblocks, nthread>>>(
                (int *)     set->cell_index_dat->data_d,
                (double *)  arg0.data_d,                // part_position,
                (double *)  arg1.data_d,                // part_weights,
                (int *)     arg2.data_d,                // part_cell_index,
                (double *)  arg3.data_d,                // part_velocity,
                (double *)  arg4.data_d,                // cell_volume,
                (double *)  arg5.data_d,                // cell_det,
                (double *)  arg6.data_d,                // cell_electric_field,
                (int *)     arg7.data_d,                // cell_connectivity_map,
                (bool *)    arg8.data_d,                // particles_injected,
                (double *)  arg9.data_d,                // dt,
                (int *)     set->particle_statuses_d,   // mark particles as MOVE_DONE, NEED_REMOVE etc...
                start, 
                end, 
                set->cells_set->size);
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());

    oppic_finalize_particle_move(set);

    cutilSafeCall(cudaFree(arg8.data_d));
    cutilSafeCall(cudaFree(arg9.data_d));
    op_mpi_set_dirtybit_cuda(nargs, args);
    cutilSafeCall(cudaDeviceSynchronize());
}

//*************************************************************************************************