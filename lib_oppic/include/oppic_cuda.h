
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

#pragma once

#include <oppic_lib.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <iostream>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/zip_iterator.h>

#define cutilSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cutilCheckMsg(msg) __cutilCheckMsg(msg, __FILE__, __LINE__)

typedef struct cudaDeviceProp cudaDeviceProp_t;

//*************************************************************************************************



//*************************************************************************************************

void __cudaSafeCall(cudaError_t err, const char *file, const int line);

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line);


void oppic_cuda_exit();

void cutilDeviceInit(int argc, char **argv);

void op_mvHostToDevice(void **map, int size);

void op_cpHostToDevice(void **data_d, void **data_h, int copy_size, int alloc_size); 

void op_upload_dat(oppic_dat dat);

void op_download_dat(oppic_dat dat);

int op_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device);

int op_mpi_halo_exchanges_cuda(oppic_set set, int nargs, oppic_arg *args);

void op_mpi_set_dirtybit_cuda(int nargs, oppic_arg *args);


void oppic_download_particle_set(oppic_set particles_set);

void oppic_upload_particle_set(oppic_set particles_set, bool realloc = false);

void oppic_create_copy_dat_to_device(oppic_dat dat);

void oppic_increase_particle_count_cuda(oppic_set particles_set);

void print_last_cuda_error();

template <class T> 
void sort_dat_according_to_index(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, int set_size)
{
    thrust::device_ptr<T> dat_dp = thrust::device_pointer_cast((T*)dat->data_d);
    thrust::device_vector<T> dat_dv(dat_dp, (dat_dp + (set_size * dat->dim)));
    thrust::device_vector<T> sorted_dat_dv(set_size * dat->dim);

    switch (dat->dim)
    {
        case 1:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(dat_dv.begin())
                ), 
                new_idx_dv.begin()), 
                set_size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(sorted_dat_dv.begin())));
            break;
        case 2:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv.begin(), 
                        (dat_dv.begin() + set_size)
                    )
                ), 
                new_idx_dv.begin()), 
                set_size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv.begin(), 
                        (sorted_dat_dv.begin() + set_size))));
            break;
        case 3:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv.begin(), 
                        (dat_dv.begin() + set_size), 
                        (dat_dv.begin() + (2 * set_size))
                    )
                ), 
                new_idx_dv.begin()), 
                set_size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv.begin(), 
                        (sorted_dat_dv.begin() + set_size), 
                        (sorted_dat_dv.begin() + (2 * set_size)))));
            break;
        case 4:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv.begin(), 
                        (dat_dv.begin() + set_size), 
                        (dat_dv.begin() + (2 * set_size)), 
                        (dat_dv.begin() + (3 * set_size))
                    )
                ), 
                new_idx_dv.begin()), 
                set_size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv.begin(), 
                        (sorted_dat_dv.begin() + set_size), 
                        (sorted_dat_dv.begin() + (2 * set_size)), 
                        (sorted_dat_dv.begin() + (3 * set_size)))));
            break;
        default:
            std::cerr << "particle_sort_cuda not implemented for dim " << dat->dim << " in dat " << dat->name << std::endl;
    }

    T* sorted_dat_dp = thrust::raw_pointer_cast(&sorted_dat_dv[0]);
    cudaMemcpy((void*)dat->data_d, (void*)sorted_dat_dp, (set_size * dat->size), cudaMemcpyDeviceToDevice);
}

void oppic_finalize_particle_move_cuda(oppic_set set);