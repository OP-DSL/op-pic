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


int dep_charge_stride_OPP_HOST_0 = -1;
int dep_charge_stride_OPP_HOST_1 = -1;

__constant__ int dep_charge_stride_OPP_CUDA_0;
__constant__ int dep_charge_stride_OPP_CUDA_1;

//user function
//*************************************************************************************************
__device__ void dep_node_charge__kernel_gpu(
    const double* part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
)
{
    atomicAdd(node_charge_den0, (part_lc[0 * dep_charge_stride_OPP_CUDA_0]));
    atomicAdd(node_charge_den1, (part_lc[1 * dep_charge_stride_OPP_CUDA_0]));
    atomicAdd(node_charge_den2, (part_lc[2 * dep_charge_stride_OPP_CUDA_0]));
    atomicAdd(node_charge_den3, (part_lc[3 * dep_charge_stride_OPP_CUDA_0]));
}

// CUDA kernel function
//*************************************************************************************************
__global__ void opp_cuda_all_DepNodeCharge(
    int *__restrict d_cell_index,
    double *__restrict dir_arg0,
    double *__restrict ind_arg1,
    const int *__restrict oppDat1Map,
    double *__restrict ind_arg2,
    double *__restrict ind_arg3,
    double *__restrict ind_arg4,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // TODO : Use shared memory to make this faster
    
    if (tid + start < end) 
    {
        int n = tid + start;

        int map0idx = d_cell_index[n];
        
        const int map1idx = oppDat1Map[map0idx + dep_charge_stride_OPP_CUDA_1 * 0];
        const int map2idx = oppDat1Map[map0idx + dep_charge_stride_OPP_CUDA_1 * 1];
        const int map3idx = oppDat1Map[map0idx + dep_charge_stride_OPP_CUDA_1 * 2];
        const int map4idx = oppDat1Map[map0idx + dep_charge_stride_OPP_CUDA_1 * 3];

// printf("n %d ci %d c_to_n %d %d %d %d\n",
//     n, map0idx, map1idx,map2idx,map3idx,map4idx);

        dep_node_charge__kernel_gpu(
            (dir_arg0 + n),
            (ind_arg1 + map1idx),
            (ind_arg2 + map2idx),
            (ind_arg3 + map3idx),
            (ind_arg4 + map4idx)
        );
    }
}

void opp_loop_all__DepositChargeOnNodes(
    opp_set set,      // particles_set
    opp_arg arg0,     // part_lc,
    opp_arg arg1,     // node_charge_den0,
    opp_arg arg2,     // node_charge_den1,
    opp_arg arg3,     // node_charge_den2,
    opp_arg arg4      // node_charge_den3,
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__DepositChargeOnNodes set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("DepCharge");

    int nargs = 5;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);

    opp_init_double_indirect_reductions_cuda(nargs, args);

    if (set_size > 0) 
    {
        dep_charge_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        dep_charge_stride_OPP_HOST_1 = args[1].size;

// opp_printf("AAAA", "dep_charge_stride_OPP_HOST_1=%d set_size=%d exec=%d nexec=%d", dep_charge_stride_OPP_HOST_1, args[1].map->from->size, args[1].map->from->exec_size, args[1].map->from->nonexec_size);
        cudaMemcpyToSymbol(dep_charge_stride_OPP_CUDA_0, &dep_charge_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(dep_charge_stride_OPP_CUDA_1, &dep_charge_stride_OPP_HOST_1, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            opp_cuda_all_DepNodeCharge <<<nblocks, nthread>>> (
                (int *)     set->mesh_relation_dat->data_d,
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,
                (int *)     args[1].map_data_d,
                (double *)  args[2].data_d,
                (double *)  args[3].data_d,
                (double *)  args[4].data_d,
                start, 
                end);
        }
    }

    opp_exchange_double_indirect_reductions_cuda(nargs, args);

    opp_complete_double_indirect_reductions_cuda(nargs, args);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("DepCharge");
}

//*************************************************************************************************