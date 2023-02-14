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

__constant__ int opDat0_MoveParticles_stride_OPPIC_CONSTANT;
__constant__ int opDat1_MoveParticles_stride_OPPIC_CONSTANT;
__constant__ int opDat2_MoveParticles_stride_OPPIC_CONSTANT;

int opDat0_MoveParticles_stride_OPPIC_HOST =- 1;
int opDat1_MoveParticles_stride_OPPIC_HOST =- 1;
int opDat2_MoveParticles_stride_OPPIC_HOST =- 1;

//user function
//*************************************************************************************************
__device__ void move_particles__kernel_gpu(
    double *pos,    
    double *vel,
    const double *ef,
    const double *dt
)
{
    double v0_add = (OP_CONST_charge / OP_CONST_mass * ef[0 * opDat2_MoveParticles_stride_OPPIC_CONSTANT] * (*dt));
    double v1_add = (OP_CONST_charge / OP_CONST_mass * ef[1 * opDat2_MoveParticles_stride_OPPIC_CONSTANT] * (*dt));
    double v2_add = (OP_CONST_charge / OP_CONST_mass * ef[2 * opDat2_MoveParticles_stride_OPPIC_CONSTANT] * (*dt));

    double v0 = (vel[0 * opDat1_MoveParticles_stride_OPPIC_CONSTANT] + v0_add);
    double v1 = (vel[1 * opDat1_MoveParticles_stride_OPPIC_CONSTANT] + v1_add);
    double v2 = (vel[2 * opDat1_MoveParticles_stride_OPPIC_CONSTANT] + v2_add);

    atomicAdd(&vel[0 * opDat1_MoveParticles_stride_OPPIC_CONSTANT], v0_add);
    atomicAdd(&vel[1 * opDat1_MoveParticles_stride_OPPIC_CONSTANT], v1_add);
    atomicAdd(&vel[2 * opDat1_MoveParticles_stride_OPPIC_CONSTANT], v2_add);
    atomicAdd(&pos[0 * opDat0_MoveParticles_stride_OPPIC_CONSTANT], (v0 * (*dt)));
    atomicAdd(&pos[1 * opDat0_MoveParticles_stride_OPPIC_CONSTANT], (v1 * (*dt)));
    atomicAdd(&pos[2 * opDat0_MoveParticles_stride_OPPIC_CONSTANT], (v2 * (*dt)));
}


// CUDA kernel function
//*************************************************************************************************
__global__ void oppic_cuda_MoveParticles(
    double *__restrict dir_arg0,
    double *__restrict dir_arg1,
    const double *__restrict dir_arg2,
    const double *__restrict dir_arg3,
    int start,
    int end,
    int set_size
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        //user-supplied kernel call
        move_particles__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (dir_arg2 + n),
            dir_arg3
        );
    }
}


//*************************************************************************************************
void oppic_par_loop_all__MoveParticles(
    oppic_set set,     // particles_set
    oppic_arg arg0,    // part_position,
    oppic_arg arg1,    // part_velocity,
    oppic_arg arg2,    // part_electric_field,
    oppic_arg arg3     // const dt 
    )
{ TRACE_ME;
    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__MoveParticles num_particles %d\n", set->size);

    int nargs = 4;
    oppic_arg args[4];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    cutilSafeCall(cudaMalloc(&(arg3.data_d), arg3.size));
    cutilSafeCall(cudaMemcpy(arg3.data_d, arg3.data, arg3.size, cudaMemcpyHostToDevice));

    int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    if (set_size > 0) 
    {
        opDat0_MoveParticles_stride_OPPIC_HOST = arg0.dat->set->size;
        opDat1_MoveParticles_stride_OPPIC_HOST = arg1.dat->set->size;
        opDat2_MoveParticles_stride_OPPIC_HOST = arg2.dat->set->size;

        cudaMemcpyToSymbol(opDat0_MoveParticles_stride_OPPIC_CONSTANT, &opDat0_MoveParticles_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat1_MoveParticles_stride_OPPIC_CONSTANT, &opDat1_MoveParticles_stride_OPPIC_HOST, sizeof(int));
        cudaMemcpyToSymbol(opDat2_MoveParticles_stride_OPPIC_CONSTANT, &opDat2_MoveParticles_stride_OPPIC_HOST, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            oppic_cuda_MoveParticles <<<nblocks, nthread>>> (
                (double *)  arg0.data_d,
                (double *)  arg1.data_d,
                (double *)  arg2.data_d,
                (double *)  arg3.data_d,
                start, 
                end, 
                set->size);
        }
    }

    cutilSafeCall(cudaFree(arg3.data_d));
    op_mpi_set_dirtybit_cuda(nargs, args);
    cutilSafeCall(cudaDeviceSynchronize());
}

//*************************************************************************************************