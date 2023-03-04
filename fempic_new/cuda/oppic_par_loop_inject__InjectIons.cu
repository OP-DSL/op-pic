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

__constant__ int injectIons_stride_OPP_CUDA_0;
__constant__ int injectIons_stride_OPP_CUDA_1;
__constant__ int injectIons_stride_OPP_CUDA_4;
__constant__ int injectIons_stride_OPP_CUDA_5;
__constant__ int injectIons_stride_OPP_CUDA_6;
__constant__ int injectIons_stride_OPP_CUDA_7;
__constant__ int injectIons_stride_OPP_CUDA_8;
__constant__ int injectIons_stride_OPP_CUDA_9;

int injectIons_stride_OPP_HOST_0 = -1;
int injectIons_stride_OPP_HOST_1 = -1;
int injectIons_stride_OPP_HOST_4 = -1;
int injectIons_stride_OPP_HOST_5 = -1;
int injectIons_stride_OPP_HOST_6 = -1;
int injectIons_stride_OPP_HOST_7 = -1;
int injectIons_stride_OPP_HOST_8 = -1;
int injectIons_stride_OPP_HOST_9 = -1;

//user function
//*************************************************************************************************
__device__ void inject_ions__kernel_gpu(
    double *part_pos,
    double *part_vel,
    int *part_cell_connectivity,
    int *cell_id, 
    double *cell_ef,
    double *iface_u,
    double *iface_v,
    double *iface_normal,
    double *node_pos,
    const double* dummy_part_random
)
{
    double a = dummy_part_random[0 * injectIons_stride_OPP_CUDA_9];
    double b = dummy_part_random[1 * injectIons_stride_OPP_CUDA_9];
    if ((a + b) > 1)  
    {
        a = (1 - a);
        b = (1 - b);
    }

    for (int i = 0; i < DIMENSIONS; i++) 
    {
        part_pos[i * injectIons_stride_OPP_CUDA_0] = a * iface_u[i * injectIons_stride_OPP_CUDA_5] + b * iface_v[i * injectIons_stride_OPP_CUDA_6] + node_pos[i * injectIons_stride_OPP_CUDA_8];

        part_vel[i * injectIons_stride_OPP_CUDA_1] = (iface_normal[i * injectIons_stride_OPP_CUDA_7] * CONST_ion_velocity_cuda);
        part_vel[i * injectIons_stride_OPP_CUDA_1] -= CONST_charge_cuda / CONST_mass_cuda * cell_ef[i * injectIons_stride_OPP_CUDA_4] * (0.5 * CONST_dt_cuda);

        // double v1 = (iface_normal[i * injectIons_stride_OPP_CUDA_7] * CONST_ion_velocity_cuda);
        // double v2 = (CONST_charge_cuda / CONST_mass_cuda);
        // double v3 = (0.5 * CONST_dt_cuda);
        // double v4 = (cell_ef[i * injectIons_stride_OPP_CUDA_4] * v3);
        // double v5 = (v2 * v4);

        // part_vel[i * injectIons_stride_OPP_CUDA_1] = (v1 - v5);
        
        // printf("%d %+2.20lE - %+2.20lE %+2.20lE\n - %+2.20lE %+2.20lE %+2.20lE %+2.20lE %+2.20lE\n\n", i,
        //     part_vel[i * injectIons_stride_OPP_CUDA_1],
        //     iface_normal[i * injectIons_stride_OPP_CUDA_7], cell_ef[i * injectIons_stride_OPP_CUDA_4],
        //     v1, 
        //     v2,
        //     v3, 
        //     v4,
        //     v5);
    }

    (*part_cell_connectivity) = (*cell_id);
}

// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_InjectIons(
    const int *__restrict mesh_relation,
    double *__restrict dir_arg0,
    double *__restrict dir_arg1,
    int *__restrict dir_arg2,
    int *__restrict ind_arg3,
    double *__restrict ind_arg4,
    const int *__restrict opDat4Map,
    double *__restrict ind_arg5,
    double *__restrict ind_arg6,
    double *__restrict ind_arg7,
    double *__restrict ind_arg8,
    double *__restrict dir_arg9,
    int start,
    int end,
    int inj_start
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const int map0idx = mesh_relation[n + inj_start]; // iface index
        const int map1idx = opDat4Map[map0idx]; // cell index

        inject_ions__kernel_gpu(
            (dir_arg0 + (n + inj_start)),
            (dir_arg1 + (n + inj_start)),
            (dir_arg2 + (n + inj_start)),
            (ind_arg3 + map0idx),
            (ind_arg4 + map1idx),
            (ind_arg5 + map0idx),
            (ind_arg6 + map0idx),
            (ind_arg7 + map0idx),
            (ind_arg8 + map0idx),
            (dir_arg9 + n)
        );
    }
}

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

    if (FP_DEBUG) printf("FEMPIC - oppic_par_loop_inject__InjectIons num_particles %d diff %d\n", set->size, set->diff);

    int nargs = 10;
    oppic_arg args[nargs] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9 };

    int set_size = oppic_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        injectIons_stride_OPP_HOST_0 = arg0.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_1 = arg1.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_4 = arg4.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_5 = arg5.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_6 = arg6.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_7 = arg7.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_8 = arg8.dat->set->set_capacity;
        injectIons_stride_OPP_HOST_9 = arg9.dat->set->set_capacity;

        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_0, &injectIons_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_1, &injectIons_stride_OPP_HOST_1, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_4, &injectIons_stride_OPP_HOST_4, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_5, &injectIons_stride_OPP_HOST_5, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_6, &injectIons_stride_OPP_HOST_6, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_7, &injectIons_stride_OPP_HOST_7, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_8, &injectIons_stride_OPP_HOST_8, sizeof(int));
        cudaMemcpyToSymbol(injectIons_stride_OPP_CUDA_9, &injectIons_stride_OPP_HOST_9, sizeof(int));

        int start     = 0;
        int end       = set->diff;
        int inj_start = (set->size - set->diff);

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_InjectIons<<<nblocks, nthread>>>(
                (int *)    set->mesh_relation_dat->data_d,
                (double *) arg0.data_d,
                (double *) arg1.data_d,
                (int *)    arg2.data_d,
                (int *)    arg3.data_d,
                (double *) arg4.data_d,
                (int *)    arg4.map_data_d,
                (double *) arg5.data_d,
                (double *) arg6.data_d,
                (double *) arg7.data_d,
                (double *) arg8.data_d,
                (double *) arg9.data_d,
                start, 
                end, 
                inj_start);
        }
    }

    oppic_mpi_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());
}