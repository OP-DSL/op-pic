
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

#include <oppic_cuda.h>

//*************************************************************************************************
__global__ void opp_cuda_AssignMeshRelation(
    int *__restrict mesh_relation,
    const int *__restrict distribution,
    int start,
    int end,
    int inj_start,
    int inlet_size
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {    
        int n = tid + start;

        for (int i = 0; i < inlet_size; i++)
        {
            if (tid < distribution[i])
            {
                // assign inlet face index as the injected particle mesh relation
                mesh_relation[n + inj_start] = i; 
                break;
            } 
        }  
    }
}


//****************************************
void opp_inc_part_count_with_distribution(opp_set particles_set, 
    int num_particles_to_insert, opp_dat iface_dist)
{
    if (OP_DEBUG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", 
        num_particles_to_insert);

    opp_profiler->start("IncPartCountWithDistribution");

    opp_dat mesh_rel_dat  = particles_set->mesh_relation_dat;

    // int nargs = 1;
    // opp_arg args[nargs];
    // args[0] = opp_get_arg(mesh_rel_dat, OP_READ);

    // int set_size = opp_mpi_halo_exchanges_grouped(particles_set, nargs, args, Device_CPU);

    // TODO : BUG What happens if the complete particle is dirty in device?

    oppic_increase_particle_count(particles_set, num_particles_to_insert);

    int nargs1 = 2;
    opp_arg args1[nargs1];

    // if iface particle distribution is dirty in device, get it to the device
    args1[0] = opp_get_arg(iface_dist, OP_READ);
    args1[1] = opp_get_arg(mesh_rel_dat, OP_WRITE);

    int set_size = opp_mpi_halo_exchanges_grouped(particles_set, nargs1, args1, Device_GPU);
    if (set_size > 0) 
    {
        int start     = 0;
        int end       = particles_set->diff;
        int inj_start = (particles_set->size - particles_set->diff);

        if (end - start > 0) 
        {
            int nthread = OPP_GPU_THREADS_PER_BLOCK;
            int nblocks = (end - start - 1) / nthread + 1;

            opp_cuda_AssignMeshRelation<<<nblocks, nthread>>>(
                (int *) mesh_rel_dat->data_d,
                (int *) iface_dist->data_d,
                start, 
                end, 
                inj_start,
                iface_dist->set->size);
        }
    }

    opp_mpi_set_dirtybit_grouped(nargs1, args1, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("IncPartCountWithDistribution");
}

