
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

#include <opp_lib.h>

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
#include <thrust/generate.h>
#include <thrust/random.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

#define cutilSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cutilCheckMsg(msg) __cutilCheckMsg(msg, __FILE__, __LINE__)

#define OPP_GPU_THREADS_PER_BLOCK 32

extern int* opp_saved_mesh_relation_d;
extern size_t opp_saved_mesh_relation_size;
extern thrust::device_vector<int> cellIdx_dv;
extern thrust::device_vector<int> i_dv;
extern char *OPP_need_remove_flags_d;

extern int *OPP_move_particle_indices_d;
extern int *OPP_move_cell_indices_d;
extern int *OPP_move_count_d;
extern thrust::device_vector<int> OPP_thrust_move_particle_indices_d;
extern thrust::device_vector<int> OPP_thrust_move_cell_indices_d;

extern int *OPP_remove_particle_indices_d;
extern int *OPP_remove_count_d;
extern thrust::device_vector<int> OPP_thrust_remove_particle_indices_d;

extern thrust::device_vector<int> ps_to_indices_dv;
extern thrust::device_vector<int> ps_from_indices_dv;
extern thrust::device_vector<int> ps_sequence_dv;

extern std::map<int, thrust::host_vector<OPP_INT>> cell_indices_hv;     // cellid in the foreign rank, arrange according to rank
extern std::map<int, thrust::host_vector<OPP_INT>> particle_indices_hv; // particle ids to send, arrange according to rank
extern std::map<int, thrust::device_vector<OPP_INT>> particle_indices_dv;
extern std::map<int, thrust::device_vector<char>> send_data;
extern std::map<int, thrust::device_vector<char>> recv_data;

// arrays for global constants and reductions
extern char *OP_reduct_h, *OP_reduct_d;

//*************************************************************************************************

void __cudaSafeCall(cudaError_t err, const char *file, const int line);

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line);


void opp_cuda_exit();

void cutilDeviceInit(int argc, char **argv);

// Copy a map from host to device
void opp_upload_map(opp_map map, bool create_new = false);

/*******************************************************************************/

void opp_halo_create();
void opp_halo_destroy();

/*******************************************************************************/

void opp_init_double_indirect_reductions_cuda(int nargs, opp_arg *args);
void opp_exchange_double_indirect_reductions_cuda(int nargs, opp_arg *args) ;
void opp_complete_double_indirect_reductions_cuda(int nargs, opp_arg *args);

/*******************************************************************************/


void print_dat_to_txtfile_mpi(opp_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name);

void print_last_cuda_error();

void opp_cpHostToDevice(void **data_d, void **data_h, size_t copy_size, size_t alloc_size = 0, 
    bool create_new = false);

void opp_create_device_arrays(opp_dat dat, bool create_new = false);

void opp_finalize_particle_move_cuda(opp_set set);

/*******************************************************************************/

void sort_dat_according_to_index_int(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
    int set_capacity, int size, bool hole_filling, int out_start_idx);
void sort_dat_according_to_index_double(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
    int set_capacity, int size, bool hole_filling, int out_start_idx);


inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_double(args, data);
#else
    (void)args;
    (void)data;
#endif
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_int(args, data);
#else
    (void)args;
    (void)data;
#endif
}

/*******************************************************************************/
// routines to resize constant/reduct arrays, if necessary

void opp_reallocReductArrays(int reduct_bytes);

void opp_mvReductArraysToDevice(int reduct_bytes);

void opp_mvReductArraysToHost(int reduct_bytes);

template <opp_access reduction, class T>
__inline__ __device__ void opp_reduction(volatile T *dat_g, T dat_l) 
{
    extern __shared__ volatile double temp2[];
    __shared__ volatile T *temp;
    temp = (T *)temp2;
    T dat_t;

    __syncthreads(); /* important to finish all previous activity */

    int tid = threadIdx.x;
    temp[tid] = dat_l;

    // first, cope with blockDim.x perhaps not being a power of 2

    __syncthreads();

    int d = 1 << (31 - __clz(((int)blockDim.x - 1)));
    // d = blockDim.x/2 rounded up to nearest power of 2

    if (tid + d < blockDim.x) {
        dat_t = temp[tid + d];

        switch (reduction) {
        case OP_INC:
        dat_l = dat_l + dat_t;
        break;
        case OP_MIN:
        if (dat_t < dat_l)
            dat_l = dat_t;
        break;
        case OP_MAX:
        if (dat_t > dat_l)
            dat_l = dat_t;
        break;
        }

        temp[tid] = dat_l;
    }

    // second, do reductions involving more than one warp

    for (d >>= 1; d > warpSize; d >>= 1) {
        __syncthreads();

        if (tid < d) {
            dat_t = temp[tid + d];

            switch (reduction) {
            case OP_INC:
                dat_l = dat_l + dat_t;
                break;
            case OP_MIN:
                if (dat_t < dat_l)
                dat_l = dat_t;
                break;
            case OP_MAX:
                if (dat_t > dat_l)
                dat_l = dat_t;
                break;
            }

            temp[tid] = dat_l;
        }
    }

    // third, do reductions involving just one warp

    __syncthreads();

    if (tid < warpSize) {
        for (; d > 0; d >>= 1) {
            __syncwarp();
            if (tid < d) {
                dat_t = temp[tid + d];

                switch (reduction) {
                case OP_INC:
                dat_l = dat_l + dat_t;
                break;
                case OP_MIN:
                if (dat_t < dat_l)
                    dat_l = dat_t;
                break;
                case OP_MAX:
                if (dat_t > dat_l)
                    dat_l = dat_t;
                break;
                }

                temp[tid] = dat_l;
            }
        }

        // finally, update global reduction variable

        if (tid == 0) {
            switch (reduction) {
            case OP_INC:
                *dat_g = *dat_g + dat_l;
                break;
            case OP_MIN:
                if (dat_l < *dat_g)
                *dat_g = dat_l;
                break;
            case OP_MAX:
                if (dat_l > *dat_g)
                *dat_g = dat_l;
                break;
            }
        }
    }
}

/*******************************************************************************/
/*
This function arranges the multi dimensional values in input array to output array according to the indices provided
    in_dat_dv - Input array with correct values
    out_dat_dv - Output array to have the sorted values
    new_idx_dv - indices of the input array to be arranged in the output array
    in_capacity - capacity of the input array (useful for multi dimensional dats)
    out_capacity - capacity of the output array (useful for multi dimensional dats)
    in_offset - start offset if the input array
    out_offset - start offset if the output array
    size - number of dat elements to arrange (usually this is the size of new_idx_dv)
    dimension - dimension of the dat
*/
template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, int in_offset, int out_offset, 
    int size, int dimension)
{
    switch (dimension)
    {
        case 1:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(in_dat_dv->begin() + in_offset)
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(out_dat_dv->begin() + out_offset)));
            break;
        case 2:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity)
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity))));
            break;
        case 3:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity), 
                        (in_dat_dv->begin() + in_offset + (2 * in_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity), 
                        (out_dat_dv->begin() + out_offset + (2 * out_capacity)))));
            break;
        case 4:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity), 
                        (in_dat_dv->begin() + in_offset + (2 * in_capacity)), 
                        (in_dat_dv->begin() + in_offset + (3 * in_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity), 
                        (out_dat_dv->begin() + out_offset + (2 * out_capacity)), 
                        (out_dat_dv->begin() + out_offset + (3 * out_capacity)))));
            break;
        default:
            std::cerr << "copy_according_to_index not implemented for dim " << dimension << std::endl;
            exit(-1);
    }
}

template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, int size, int dimension)
{
    copy_according_to_index<T>(in_dat_dv, out_dat_dv, new_idx_dv, in_capacity, out_capacity, 
        0, 0, size, dimension);
}

/*******************************************************************************/