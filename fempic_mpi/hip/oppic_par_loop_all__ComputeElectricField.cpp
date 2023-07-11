#include "hip/hip_runtime.h"
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


__constant__ int computeElectricField_stride_OPP_CUDA_0;
__constant__ int computeElectricField_stride_OPP_CUDA_1;

int computeElectricField_stride_OPP_HOST_0 = -1;
int computeElectricField_stride_OPP_HOST_1 = -1;

//user function
//*************************************************************************************************
__device__ void compute_electric_field__kernel_gpu(
    double *cell_electric_field,             
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
)
{
    for (int dim = 0; dim < DIMENSIONS; dim++)
    {
        cell_electric_field[dim * computeElectricField_stride_OPP_CUDA_0] -= (cell_shape_deriv[(0 * DIMENSIONS + dim) * computeElectricField_stride_OPP_CUDA_1] * (*node_potential0));
        cell_electric_field[dim * computeElectricField_stride_OPP_CUDA_0] -= (cell_shape_deriv[(1 * DIMENSIONS + dim) * computeElectricField_stride_OPP_CUDA_1] * (*node_potential1));
        cell_electric_field[dim * computeElectricField_stride_OPP_CUDA_0] -= (cell_shape_deriv[(2 * DIMENSIONS + dim) * computeElectricField_stride_OPP_CUDA_1] * (*node_potential2));
        cell_electric_field[dim * computeElectricField_stride_OPP_CUDA_0] -= (cell_shape_deriv[(3 * DIMENSIONS + dim) * computeElectricField_stride_OPP_CUDA_1] * (*node_potential3));
    }    
}

// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_ComputeElectricField(
    double *__restrict dir_arg0,
    const double *__restrict dir_arg1,
    const double *__restrict ind_arg2,
    const int *__restrict opDat2Map,
    const double *__restrict ind_arg3,
    const double *__restrict ind_arg4,
    const double *__restrict ind_arg5,
    int start,
    int end,
    int set_size) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const int map1idx = opDat2Map[n + set_size * 0];
        const int map2idx = opDat2Map[n + set_size * 1];
        const int map3idx = opDat2Map[n + set_size * 2];
        const int map4idx = opDat2Map[n + set_size * 3];

        //user-supplied kernel call
        compute_electric_field__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (ind_arg2 + map1idx),
            (ind_arg2 + map2idx),
            (ind_arg2 + map3idx),
            (ind_arg2 + map4idx)
        );
    }
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
{ 

    if (FP_DEBUG) printf("FEMPIC - oppic_par_loop_all__ComputeElectricField set_size %d\n", set->size);

    const int nargs = 6;
    oppic_arg args[nargs] = { arg0, arg1, arg2, arg3, arg4, arg5 };

    int set_size = oppic_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        computeElectricField_stride_OPP_HOST_0 = arg0.dat->set->set_capacity;
        computeElectricField_stride_OPP_HOST_1 = arg1.dat->set->set_capacity;

        hipMemcpyToSymbol(HIP_SYMBOL(computeElectricField_stride_OPP_CUDA_0), &computeElectricField_stride_OPP_HOST_0, sizeof(int));
        hipMemcpyToSymbol(HIP_SYMBOL(computeElectricField_stride_OPP_CUDA_1), &computeElectricField_stride_OPP_HOST_1, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            // int nthread = GPU_THREADS_PER_BLOCK;
            int nthread = opp_params->get<INT>("opp_threads_per_block");
            int nblocks = (end - start - 1) / nthread + 1;

            // oppic_cuda_ComputeElectricField <<<nblocks, nthread>>> (
            //     (double *)  arg0.data_d,
            //     (double *)  arg1.data_d,       
            //     (double *)  arg2.data_d,
            //     (int *)     arg2.map_data_d,
            //     (double *)  arg3.data_d,
            //     (double *)  arg4.data_d,
            //     (double *)  arg5.data_d,
            //     start, 
            //     end,
            //     arg2.map->from->size);
            
            hipLaunchKernelGGL(oppic_cuda_ComputeElectricField, 
                dim3(nblocks),
                dim3(nthread),
                0, 0,
                (double *)  arg0.data_d,
                (double *)  arg1.data_d,       
                (double *)  arg2.data_d,
                (int *)     arg2.map_data_d,
                (double *)  arg3.data_d,
                (double *)  arg4.data_d,
                (double *)  arg5.data_d,
                start, 
                end,
                arg2.map->from->size);
        }
    }

    oppic_mpi_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());
}

//*************************************************************************************************
