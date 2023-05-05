
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

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
// #include <device_launch_parameters.h>

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

//*************************************************************************************************

void __cudaSafeCall(hipError_t err, const char *file, const int line);

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line);


void oppic_cuda_exit();

void cutilDeviceInit(int argc, char **argv);

void oppic_upload_dat(oppic_dat dat);

void oppic_download_dat(oppic_dat dat);

void oppic_download_particle_set(oppic_set particles_set);

void oppic_upload_particle_set(oppic_set particles_set, bool realloc = false);

int oppic_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device);

void oppic_mpi_set_dirtybit_grouped(int nargs, oppic_arg *args, DeviceType device);

void print_last_cuda_error();

void oppic_cpHostToDevice(void **data_d, void **data_h, int copy_size, int alloc_size = 0, bool create_new = false);

void oppic_create_device_arrays(oppic_dat dat, bool create_new = false);

void oppic_finalize_particle_move_cuda(oppic_set set);

void sort_dat_according_to_index_int(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, int set_capacity, int size);
void sort_dat_according_to_index_double(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, int set_capacity, int size);

template <class T> 
void copy_according_to_index(thrust::device_vector<T>* dat_dv, thrust::device_vector<T>* sorted_dat_dv, const thrust::device_vector<int>& new_idx_dv, int set_capacity, int size, int dimension)
{
    switch (dimension)
    {
        case 1:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(dat_dv->begin())
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(sorted_dat_dv->begin())));
            break;
        case 2:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv->begin(), 
                        (dat_dv->begin() + set_capacity)
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv->begin(), 
                        (sorted_dat_dv->begin() + set_capacity))));
            break;
        case 3:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv->begin(), 
                        (dat_dv->begin() + set_capacity), 
                        (dat_dv->begin() + (2 * set_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv->begin(), 
                        (sorted_dat_dv->begin() + set_capacity), 
                        (sorted_dat_dv->begin() + (2 * set_capacity)))));
            break;
        case 4:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        dat_dv->begin(), 
                        (dat_dv->begin() + set_capacity), 
                        (dat_dv->begin() + (2 * set_capacity)), 
                        (dat_dv->begin() + (3 * set_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        sorted_dat_dv->begin(), 
                        (sorted_dat_dv->begin() + set_capacity), 
                        (sorted_dat_dv->begin() + (2 * set_capacity)), 
                        (sorted_dat_dv->begin() + (3 * set_capacity)))));
            break;
        default:
            std::cerr << "particle_sort_cuda not implemented for dim " << dimension << std::endl;
            exit(-1);
    }
}
