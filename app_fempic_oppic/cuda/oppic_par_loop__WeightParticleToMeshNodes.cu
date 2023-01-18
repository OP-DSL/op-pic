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

__constant__ int opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT;

int opDat0_WeightParticleToMeshNodes_stride_OPPIC_HOST =- 1;


//user function
//*************************************************************************************************
__device__ void weight_particle_to_mesh_nodes__kernel_gpu(
    const double* part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3,
    const double *node_volume0,
    const double *node_volume1,
    const double *node_volume2,
    const double *node_volume3
)
{
    (*node_charge_den0) += (part_lc[0 * opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT] * (OP_CONST_CUDA_spwt / (*node_volume0)));
    (*node_charge_den1) += (part_lc[1 * opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT] * (OP_CONST_CUDA_spwt / (*node_volume1)));
    (*node_charge_den2) += (part_lc[2 * opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT] * (OP_CONST_CUDA_spwt / (*node_volume2)));
    (*node_charge_den3) += (part_lc[3 * opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT] * (OP_CONST_CUDA_spwt / (*node_volume3)));
}


// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_WeightParticleToMeshNodes(
    const int *__restrict d_cell_index,
    const double *__restrict dir_arg0,
    double *__restrict ind_arg1,
    const int *__restrict opDat1Map,
    const double *__restrict ind_arg5,
    const int *__restrict opDat5Map,
    int start,
    int end,
    int set_size
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;
        
        int map0idx = d_cell_index[n];

        int map1idx = opDat1Map[map0idx + set_size * 0];
        int map2idx = opDat1Map[map0idx + set_size * 1];
        int map3idx = opDat1Map[map0idx + set_size * 2];
        int map4idx = opDat1Map[map0idx + set_size * 3];

        int map5idx = opDat5Map[map0idx + set_size * 0];
        int map6idx = opDat5Map[map0idx + set_size * 1];
        int map7idx = opDat5Map[map0idx + set_size * 2];
        int map8idx = opDat5Map[map0idx + set_size * 3];

        double arg1_l = ZERO_double;
        double arg2_l = ZERO_double;
        double arg3_l = ZERO_double;
        double arg4_l = ZERO_double;

        //user-supplied kernel call
        weight_particle_to_mesh_nodes__kernel_gpu(
            (dir_arg0 + n),
            &(arg1_l),
            &(arg2_l),
            &(arg3_l),
            &(arg4_l),
            (ind_arg5 + map5idx),
            (ind_arg5 + map6idx),
            (ind_arg5 + map7idx),
            (ind_arg5 + map8idx)
        );

        atomicAdd(&(ind_arg1[map1idx]), arg1_l);
        atomicAdd(&(ind_arg1[map2idx]), arg2_l);
        atomicAdd(&(ind_arg1[map3idx]), arg3_l);
        atomicAdd(&(ind_arg1[map4idx]), arg4_l);
    }
}

//*************************************************************************************************
void oppic_par_loop_all__WeightParticleToMeshNodes(
    oppic_set set,         // particles_set
    oppic_arg arg0,        // particle_lc
    oppic_arg arg1,        // node_charge_density
    oppic_arg arg2,        // node_charge_density
    oppic_arg arg3,        // node_charge_density
    oppic_arg arg4,        // node_charge_density
    oppic_arg arg5,        // node_volumes
    oppic_arg arg6,        // node_volumes
    oppic_arg arg7,        // node_volumes
    oppic_arg arg8         // node_volumes        
    )
{ TRACE_ME;
    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__WeightParticleToMeshNodes num_particles %d\n", set->size);

    int nargs = 9;
    oppic_arg args[9];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;
    args[7] = arg7;
    args[8] = arg8;

    int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        opDat0_WeightParticleToMeshNodes_stride_OPPIC_HOST = arg0.dat->set->size;
        cudaMemcpyToSymbol(opDat0_WeightParticleToMeshNodes_stride_OPPIC_CONSTANT, &opDat0_WeightParticleToMeshNodes_stride_OPPIC_HOST, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_WeightParticleToMeshNodes <<<nblocks, nthread>>> (
                (int *)     set->cell_index_dat->data_d,
                (double *)  arg0.data_d,
                (double *)  arg1.data_d,
                (int *)     arg1.map_data_d,
                (double *)  arg5.data_d,
                (int *)     arg5.map_data_d,
                start, 
                end, 
                arg1.map->from->size);
        }
    }

    op_mpi_set_dirtybit_cuda(nargs, args);
    cutilSafeCall(cudaDeviceSynchronize());
}

//*************************************************************************************************