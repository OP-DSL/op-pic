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


inline void calculate_injection_distribution(
    int* injected_total,
    double* face_area,
    int* particle_distribution,
    double* remainder 
)
{   
    /*number of real ions per sec, given prescribed density and velocity*/
    double num_per_sec = CONST_plasma_den * CONST_ion_velocity * (*face_area);

    /*number of ions to generate in this time step*/
    double num_real = num_per_sec * CONST_dt;

    /*fraction number of macroparticles*/
    double fnum_mp = num_real / CONST_spwt + (*remainder);

    /*integer number of macroparticles*/
    int num_mp = (int)fnum_mp;

    /*update reminder*/
    (*remainder) = fnum_mp - num_mp;

    (*injected_total) += num_mp;

    (*particle_distribution) = (*injected_total);
}

#ifdef USE_CPU_ASSIGN

void oppic_inject__Increase_particle_count(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
)
{ TRACE_ME;

    if (FP_DEBUG) printf("FEMPIC - oppic_inject__Increase_particle_count set_size %d diff %d\n", set->size, set->diff);

    int nargs = 5; // Add one more to have mesh_relation arg
    oppic_arg args[nargs];
    args[0] = arg0; 
    args[1] = arg1; 
    args[2] = arg2; 
    args[3] = arg3; 
    args[4] = oppic_arg_dat(particles_set->mesh_relation_dat, OP_RW);
    
    particles_set->mesh_relation_dat->dirty_hd = Dirty::Host; // make mesh relation dirty and download new data from device

    int set_size = oppic_mpi_halo_exchanges_grouped(set, nargs, args, Device_CPU);

    for (int i = 0; i < set_size; i++)
    {   
        calculate_injection_distribution(
            ((int *)arg0.data),
            &((double *)arg1.data)[i],
            &((int *)arg2.data)[i],
            ((double *)arg3.data) 
        );
    }

    oppic_increase_particle_count(particles_set, *((int *)arg0.data));

    oppic_dat mesh_rel_dat  = particles_set->mesh_relation_dat;
    int* part_mesh_relation = (int *)mesh_rel_dat->data;
    int* distribution       = (int *)arg2.data;

    int start = (particles_set->size - particles_set->diff);
    int j = 0;

    for (int i = 0; i < particles_set->diff; i++)
    {
        if (i >= distribution[j]) j++; // check whether it is j or j-1    
        part_mesh_relation[start + i] = j;
    }  

    oppic_mpi_set_dirtybit_grouped(nargs, args, Device_CPU);

    int size = start * mesh_rel_dat->size;
    char* inj_start_data   = (mesh_rel_dat->data + size);
    char* inj_start_data_d = (mesh_rel_dat->data_d + size);

    oppic_cpHostToDevice((void **)&(inj_start_data_d), (void **)&(inj_start_data), (mesh_rel_dat->size * particles_set->diff));

    mesh_rel_dat->dirty_hd = Dirty::NotDirty;
}

#else

//*************************************************************************************************
__global__ void oppic_cuda_AssignMeshRelation(
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
                mesh_relation[n + inj_start] = i;
                break;
            } 
        }  
    }
}


void oppic_inject__Increase_particle_count(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
)
{ TRACE_ME;

    if (FP_DEBUG) printf("FEMPIC - oppic_inject__Increase_particle_count set_size %d diff %d\n", set->size, set->diff);

    const int nargs = 4; // Add one more to have mesh_relation arg
    oppic_arg args[nargs] = { arg0, arg1, arg2, arg3 };

    int set_size = oppic_mpi_halo_exchanges_grouped(set, nargs, args, Device_CPU);

    for (int i = 0; i < set_size; i++)
    {   
        calculate_injection_distribution(
            ((int *)arg0.data),
            &((double *)arg1.data)[i],
            &((int *)arg2.data)[i],
            ((double *)arg3.data) 
        );
    }

    oppic_increase_particle_count(particles_set, *((int *)arg0.data));

    oppic_mpi_set_dirtybit_grouped(nargs, args, Device_CPU);

    arg2.acc = OP_READ; // if not, iface_inj_part_dist will not get uploaded to device
    oppic_dat mesh_rel_dat  = particles_set->mesh_relation_dat;

    const int nargs1 = 2; // Add one more to have mesh_relation arg
    oppic_arg args1[nargs1] = { arg2, oppic_arg_dat(mesh_rel_dat, OP_WRITE) };

    set_size = oppic_mpi_halo_exchanges_grouped(particles_set, nargs1, args1, Device_GPU);
    if (set_size > 0) 
    {
        int start     = 0;
        int end       = particles_set->diff;
        int inj_start = (particles_set->size - particles_set->diff);

        if (end - start > 0) 
        {
            // int nthread = GPU_THREADS_PER_BLOCK;
            int nthread = opp_params->get<INT>("opp_threads_per_block");
            int nblocks = (end - start - 1) / nthread + 1;

            // oppic_cuda_AssignMeshRelation<<<nblocks, nthread>>>(
            //     (int *) mesh_rel_dat->data_d,
            //     (int *) arg2.data_d,
            //     start, 
            //     end, 
            //     inj_start,
            //     set->size);


            hipLaunchKernelGGL(oppic_cuda_AssignMeshRelation, 
                dim3(nblocks),
                dim3(nthread),
                0, 0,
                (int *) mesh_rel_dat->data_d,
                (int *) arg2.data_d,
                start, 
                end, 
                inj_start,
                set->size);
        }
    }

    oppic_mpi_set_dirtybit_grouped(nargs1, args1, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());
}

#endif




