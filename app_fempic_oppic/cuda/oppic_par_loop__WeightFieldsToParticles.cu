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

__constant__ int opDat0_WeightFieldsToParticles_stride_OPPIC_CONSTANT;
__constant__ int opDat1_WeightFieldsToParticles_stride_OPPIC_CONSTANT;

int opDat0_WeightFieldsToParticles_stride_OPPIC_HOST =- 1;
int opDat1_WeightFieldsToParticles_stride_OPPIC_HOST =- 1;


//user function
//*************************************************************************************************
__device__ void weight_fields_to_particles__kernel_gpu(
    double *part_ef,
    const double *cell_ef
)
{
    part_ef[0 * opDat0_WeightFieldsToParticles_stride_OPPIC_CONSTANT] = cell_ef[0 * opDat1_WeightFieldsToParticles_stride_OPPIC_CONSTANT];
    part_ef[1 * opDat0_WeightFieldsToParticles_stride_OPPIC_CONSTANT] = cell_ef[1 * opDat1_WeightFieldsToParticles_stride_OPPIC_CONSTANT];
    part_ef[2 * opDat0_WeightFieldsToParticles_stride_OPPIC_CONSTANT] = cell_ef[2 * opDat1_WeightFieldsToParticles_stride_OPPIC_CONSTANT];
}


// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_WeightFieldsToParticles(
    const int *__restrict d_cell_index,
    double *__restrict dir_arg0,
    const double *__restrict ind_arg1,
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

        //user-supplied kernel call
        weight_fields_to_particles__kernel_gpu(
            (dir_arg0 + n),
            (ind_arg1 + map0idx)
        );
    }
}


//*************************************************************************************************
void oppic_par_loop_all__WeightFieldsToParticles(
    oppic_set set,     // particles_set
    oppic_arg arg0,    // particle_ef
    oppic_arg arg1     // cell_electric_field
    )
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__WeightFieldsToParticles num_particles %d\n", set->size);

    int nargs = 2;
    oppic_arg args[2];

    args[0] = arg0;
    args[1] = arg1;

    int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        opDat0_WeightFieldsToParticles_stride_OPPIC_HOST = arg0.dat->set->size;
        opDat1_WeightFieldsToParticles_stride_OPPIC_HOST = arg1.dat->set->size;
        
        cudaMemcpyToSymbol(opDat0_WeightFieldsToParticles_stride_OPPIC_CONSTANT, &opDat0_WeightFieldsToParticles_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat1_WeightFieldsToParticles_stride_OPPIC_CONSTANT, &opDat1_WeightFieldsToParticles_stride_OPPIC_HOST, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_WeightFieldsToParticles <<<nblocks, nthread>>> (
                (int *)     set->cell_index_dat->data_d,
                (double *)  arg0.data_d,
                (double *)  arg1.data_d,
                start, 
                end, 
                set->size);
        }
    }

    op_mpi_set_dirtybit_cuda(nargs, args);
    cutilSafeCall(cudaDeviceSynchronize());
}

//*************************************************************************************************