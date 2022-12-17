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


//user function
//*************************************************************************************************
__device__ void reset_ion_density__kernel_gpu(
    double *ion_den
)
{
    *ion_den = ZERO_double;
}


// CUDA kernel function
//*************************************************************************************************
__global__ void op_cuda_ResetIonDensity(
    double *__restrict dir_arg0,
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
        reset_ion_density__kernel_gpu(
            (dir_arg0 + n)
        );
    }
}


//*************************************************************************************************
void op_par_loop_all__ResetIonDensity(
    op_set set,     // nodes_set
    op_arg arg0     // node_charge_density
    )
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - op_par_loop_all__ResetIonDensity num_nodes %d\n", set->size);

    int start = 0;
    int end   = set->size;

    if (end - start > 0) 
    {
        int nthread = GPU_THREADS_PER_BLOCK;
        int nblocks = (end - start - 1) / nthread + 1;

        op_cuda_ResetIonDensity <<<nblocks, nthread>>> (
            (double *)  arg0.data_d,
            start, 
            end, 
            set->size);
    }  
}

//*************************************************************************************************