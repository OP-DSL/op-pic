
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

#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <dpct/dpl_utils.hpp>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/random>

#include <opp_lib.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

constexpr bool debug_mem = false;

#define cutilSafeCall(err) // these need to be remmoved
#define cutilCheckMsg(msg) // these need to be remmoved

#define OPP_GPU_THREADS_PER_BLOCK 32

#define OPP_PARTICLE_MOVE_DONE { m.move_status = OPP_MOVE_DONE; }
#define OPP_PARTICLE_NEED_MOVE { m.move_status = OPP_NEED_MOVE; }
#define OPP_PARTICLE_NEED_REMOVE { m.move_status = OPP_NEED_REMOVE; }
#define OPP_DO_ONCE (m.iteration_one)
#define OPP_MOVE_RESET_FLAGS { m.move_status = OPP_MOVE_DONE; m.iteration_one = true; }

extern int* opp_saved_mesh_relation_d;
extern size_t opp_saved_mesh_relation_size;

extern dpct::device_vector<int> ps_cell_index_dv;
extern dpct::device_vector<int> ps_swap_indices_dv;
// extern char *OPP_need_remove_flags_d;

extern int *OPP_move_particle_indices_d;
extern int *OPP_move_cell_indices_d;
extern int *OPP_move_count_d;
// extern dpct::device_vector<int> OPP_thrust_move_particle_indices_d;
// extern dpct::device_vector<int> OPP_thrust_move_cell_indices_d;

extern int *OPP_remove_particle_indices_d;
// extern int *OPP_remove_count_d;
// extern dpct::device_vector<int> OPP_thrust_remove_particle_indices_d;

// extern dpct::device_vector<int> ps_to_indices_dv;
extern dpct::device_vector<int> hf_from_indices_dv;
extern dpct::device_vector<int> hf_sequence_dv;

extern std::map<int, std::vector<OPP_INT>>
    cell_indices_hv; // cellid in the foreign rank, arrange according to rank
extern std::map<int, std::vector<OPP_INT>>
    particle_indices_hv; // particle ids to send, arrange according to rank
extern std::map<int, dpct::device_vector<OPP_INT>> particle_indices_dv;
extern std::map<int, dpct::device_vector<char>> send_data;
extern std::map<int, dpct::device_vector<char>> recv_data;

// arrays for global constants and reductions
extern int OPP_consts_bytes;
extern int OPP_reduct_bytes;
extern char *OPP_reduct_h, *OPP_reduct_d;
extern char *OPP_consts_h, *OPP_consts_d;

extern std::vector<char*> opp_consts;

//*************************************************************************************************

void __cudaSafeCall(dpct::err0 err, const char *file, const int line);

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line);

void opp_sycl_init(int argc, char **argv);
void opp_sycl_exit();

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

void opp_copy_host_to_device(void **data_d, void **data_h, size_t copy_size, 
    size_t alloc_size = 0, bool create_new = false);

void opp_create_dat_device_arrays(opp_dat dat, bool create_new = false);

void opp_finalize_particle_move_cuda(opp_set set);

void particle_sort_device(opp_set set, bool hole_filling);

/*******************************************************************************/

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

void opp_reallocConstArrays(int consts_bytes);

void opp_mvConstArraysToDevice(int consts_bytes);

void opp_mvConstArraysToHost(int consts_bytes);

template <opp_access reduction, int intel, class T, class out_acc, class local_acc>
void opp_reduction(out_acc dat_g, int offset, T dat_l, local_acc temp, sycl::nd_item<1> &item) {
    T dat_t;

    /* important to finish all previous activity */
    item.barrier(sycl::access::fence_space::local_space); 

    size_t tid = item.get_local_id(0);
    temp[tid] = dat_l;

    for (size_t d = item.get_local_range(0) / 2; d > 0; d >>= 1) {
        item.barrier(sycl::access::fence_space::local_space);
        if (tid < d) {
        dat_t = temp[tid + d];

        switch (reduction) {
        case OPP_INC:
            dat_l = dat_l + dat_t;
            break;
        case OPP_MIN:
            if (dat_t < dat_l)
                dat_l = dat_t;
            break;
        case OPP_MAX:
            if (dat_t > dat_l)
                dat_l = dat_t;
            break;
        }
        temp[tid] = dat_l;
        }
    }

    if (tid == 0) {
        switch (reduction) {
        case OPP_INC:
            dat_g[offset] = dat_g[offset] + dat_l;
            break;
        case OPP_MIN:
            if (dat_l < dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        case OPP_MAX:
            if (dat_l > dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        }
    }
}

/*******************************************************************************/
/*
This function arranges the multi dimensional values in input array to output array 
according to the indices provided
    in_dat_dv - Input array with correct values
    out_dat_dv - Output array to have the sorted values
    new_idx_dv - indices of the input array to be arranged in the output array
    in_capacity - capacity of the input array (useful for multi dimensional dats)
    out_capacity - capacity of the output array (useful for multi dimensional dats)
    in_offset - start offset if the input array
    out_offset - start offset if the output array
    size - number of dat elements to arrange (usually this is the size of new_idx_dv)
    dim - dimension of the dat
*/
template <class T>
void copy_according_to_index(const T *in_dat_dv, T *out_dat_dv, const OPP_INT* new_idx_dv,
        int in_capacity, int out_capacity, int in_offset, int out_offset, int size, int dim)
{
    switch (dim)
    {
        case 1:
            std::copy_n(oneapi::dpl::execution::make_device_policy(*opp_queue),
                oneapi::dpl::make_permutation_iterator(
                    oneapi::dpl::make_zip_iterator(std::make_tuple(
                        in_dat_dv + in_offset)),
                    new_idx_dv),
                size,
                oneapi::dpl::make_zip_iterator(std::make_tuple(
                    out_dat_dv + out_offset)
                ));
            break;
        case 2:
            std::copy_n(oneapi::dpl::execution::make_device_policy(*opp_queue),
                oneapi::dpl::make_permutation_iterator(
                    oneapi::dpl::make_zip_iterator(std::make_tuple(
                        (in_dat_dv + in_offset),
                        (in_dat_dv + in_offset + in_capacity))),
                    new_idx_dv),
                size,
                oneapi::dpl::make_zip_iterator(std::make_tuple(
                    (out_dat_dv + out_offset),
                    (out_dat_dv + out_offset + out_capacity))));
            break;
        case 3:
            std::copy_n(oneapi::dpl::execution::make_device_policy(*opp_queue),
                oneapi::dpl::make_permutation_iterator(
                    oneapi::dpl::make_zip_iterator(std::make_tuple(
                        (in_dat_dv + in_offset),
                        (in_dat_dv + in_offset + in_capacity))),
                    new_idx_dv),
                size,
                oneapi::dpl::make_zip_iterator(std::make_tuple(
                    (out_dat_dv + out_offset),
                    (out_dat_dv + out_offset + out_capacity))));
            break;
        case 4:
            std::copy_n(oneapi::dpl::execution::make_device_policy(*opp_queue),
                oneapi::dpl::make_permutation_iterator(
                    oneapi::dpl::make_zip_iterator(std::make_tuple(
                        (in_dat_dv + in_offset),
                        (in_dat_dv + in_offset + in_capacity))),
                    new_idx_dv),
                size,
                oneapi::dpl::make_zip_iterator(std::make_tuple(
                    (out_dat_dv + out_offset),
                    (out_dat_dv + out_offset + out_capacity))));
            break;
        default:
            std::cerr << "copy_according_to_index not implemented for dim " << dim << std::endl;
            exit(-1);
    }
}

template <class T>
void copy_according_to_index(dpct::device_vector<T> *in_dat_dv,
                             dpct::device_vector<T> *out_dat_dv,
                             const dpct::device_vector<int> &new_idx_dv,
                             int in_capacity, int out_capacity, int in_offset,
                             int out_offset, int size, int dim)
{
    copy_according_to_index<T>(
                        dpct::get_raw_pointer(in_dat_dv->data()), 
                        dpct::get_raw_pointer(out_dat_dv->data()), 
                        dpct::get_raw_pointer(new_idx_dv.data()),
                        in_capacity, out_capacity, in_offset, out_offset, size, dim);
}

template <class T>
void copy_according_to_index(const T *in_dat_dv, T *out_dat_dv, const OPP_INT* new_idx_dv,
        int in_capacity, int out_capacity, int size, int dim)
{
    copy_according_to_index(in_dat_dv, out_dat_dv, new_idx_dv, in_capacity, out_capacity, 
        0, 0, size, dim);
}

template <class T>
void copy_according_to_index(dpct::device_vector<T> *in_dat_dv,
                             dpct::device_vector<T> *out_dat_dv,
                             const dpct::device_vector<int> &new_idx_dv,
                             int in_capacity, int out_capacity, int size,
                             int dim)
{
    copy_according_to_index<T>(in_dat_dv, out_dat_dv, new_idx_dv, in_capacity, out_capacity, 
        0, 0, size, dim);
}

// template <class T>
// void opp_device_memset(T* data, T value, size_t size) {
//     opp_queue->submit([&](sycl::handler& cgh) {
//         cgh.parallel_for<class memset_kernel>(sycl::range<1>(size), [=](sycl::id<1> idx) {
//             data[idx] = value;
//         });
//     }).wait(); // Wait for the kernel to finish
// }

/*******************************************************************************/
class opp_mem {

public:

    // Allocate host memory
    template <typename T>
    static T* host_malloc(size_t count) {
        return (T*)malloc(count * sizeof(T));
    }

    // Free host memory
    template <typename T>
    static void host_free(T* ptr) {
        if (ptr)
            free(ptr);
        ptr = nullptr;
    }

    // Reallocate host memory
    template <typename T>
    static void host_realloc(T*& ptr, size_t new_size) {
        T* tmp_ptr = (T*)realloc(ptr, new_size);
        ptr = tmp_ptr;
    }

    // Allocate device memory
    template <typename T>
    static T* dev_malloc(size_t count, sycl::queue* q = opp_queue) {
        if (count <= 0) return nullptr;
        try {
            T* ptr = sycl::malloc_device<T>(count * sizeof(T), *q);
            if (debug_mem) opp_printf("dev_malloc", "[%p][%zu]", ptr, count);
            return ptr;
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("dev_malloc: ") + e.what());
        }
    }

    // Free device memory
    template <typename T>
    static void dev_free(T*& ptr, sycl::queue* q = opp_queue) {
        if (ptr) {
            if (debug_mem) opp_printf("dev_free", "[%p]", ptr);
            sycl::free(ptr, *q);
        }
    }

    // Copy memory from one device pointer to another
    template <typename T>
    static void dev_memcpy(T*& dst, const T* src, size_t cpy_count, sycl::queue* q = opp_queue) {
        if (debug_mem) opp_printf("dev_memcpy", "[%p]->[%p] cpy_count[%zu]", src, dst, cpy_count);
        q->memcpy(dst, src, cpy_count * sizeof(T)).wait();
    }

    // Resize device memory
    template <typename T>
    static void dev_realloc(T*& ptr, size_t& current_size, const size_t& new_size, sycl::queue* q = opp_queue) {
        if (new_size <= 0) 
            throw std::runtime_error("dev_realloc: New Realloc size invalid - " + std::to_string(new_size));
        T* new_ptr = opp_mem::dev_malloc<T>(new_size, q);
        if (debug_mem) opp_printf("dev_realloc", "created [%p] old [%p]", new_ptr, ptr);
        if (ptr) {
            const size_t copy_size = std::min(current_size, new_size);
            opp_mem::dev_memcpy<T>(new_ptr, ptr, copy_size, q);
            opp_mem::dev_free<T>(ptr, q);
            current_size = new_size;
            if (debug_mem) opp_printf("dev_realloc", "[%p]->[%p] cpy_count[%zu]", ptr, new_ptr, copy_size);
        }
        ptr = new_ptr;
    }

    // Resize device memory (only increasing size)
    template <typename T>
    static void dev_resize(T*& ptr, size_t& current_size, const size_t& new_size, sycl::queue* q = opp_queue) {
        if (debug_mem) opp_printf("dev_resize", "[%p] %zu -> %zu", ptr, current_size, new_size);
        if (new_size > current_size) {
            opp_mem::dev_realloc<T>(ptr, current_size, new_size, q);
        }
    }

    // initialize device memory with a specific value
    template <typename T>
    static void dev_memset(T* ptr, size_t count, T value, sycl::queue* q = opp_queue) {
        q->fill(ptr, value, count).wait();
    }

    // Allocate and initialize device memory with a specific value
    template <typename T>
    static T* dev_malloc_set(size_t count, T value, sycl::queue* q = opp_queue) {
        T* ptr = opp_mem::dev_malloc<T>(count, q);
        opp_mem::dev_memset<T>(ptr, count, value, q);
        return ptr;
    }

    // Copy data from host to device, create new device arrays if requested
    template <typename T>
    static void copy_host_to_dev(T*& data_d, const T *data_h, size_t copy_count, sycl::queue* q = opp_queue, 
                        bool no_wait = false, bool create_new = false, size_t alloc_count = 0) {
        try {
            if (create_new) {
                if (data_d != nullptr)  
                    opp_mem::dev_free<T>(data_d, q);
                data_d = opp_mem::dev_malloc<T>(alloc_count, q);
            }
            q->memcpy(data_d, data_h, copy_count * sizeof(T));
            if (!no_wait) 
                opp_queue->wait(); 
            if (debug_mem) opp_printf("copy_host_to_dev", "[%p]->[%p] copy_count[%zu]", data_h, data_d, copy_count);
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("opp_mem::copy_host_to_dev: ") + e.what());
        }
    }

    // Copy data from device to host, no dot create new host arrays since it can be allocated differently, 
    // like malloc, new, stack, std::vector<>, hence free mechanism is unknown
    template <typename T>
    static void copy_dev_to_host(T* data_h, const T *data_d, size_t copy_count, sycl::queue* q = opp_queue, 
                        bool no_wait = false) {
        try {
            q->memcpy(data_h, data_d, copy_count * sizeof(T));
            if (!no_wait) 
                opp_queue->wait();
            if (debug_mem) opp_printf("copy_dev_to_host", "[%p]->[%p] copy_count[%zu]", data_d, data_h, copy_count);
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("opp_mem::copy_dev_to_host: ") + e.what());
        }
    }
};

//****************************************
template <typename T>
void copy_from(
    const T* in_dat_d, T* out_dat_d, 
    const int* from_idx_map, 
    const int in_stride, const int out_stride, 
    const int in_offset, const int out_offset, 
    const int dim, const int size,
    const sycl::nd_item<1> &item) 
{
    const int tid = item.get_global_linear_id();
    if (tid < size) {
        const int idx = from_idx_map[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + tid + d * in_stride] = 
                            in_dat_d[in_offset + idx + d * out_stride];
        }
    }
}

//****************************************
template <typename T>
void copy_from_to(
    const T* in_dat_d, T* out_dat_d, 
    const int* from_idx_map, const int* to_idx_map, 
    const int in_stride, const int out_stride, 
    const int in_offset, const int out_offset, 
    const int dim, const int size,
    const sycl::nd_item<1> &item) 
{
    const int tid = item.get_global_linear_id();
    if (tid < size) {
        const int f_idx = from_idx_map[tid];
        const int t_idx = to_idx_map[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + t_idx + d * in_stride] = 
                            in_dat_d[in_offset + f_idx + d * out_stride];
        }
    }
}

template <class T>
void opp_register_const(T*& ptr, const size_t count) {
    if (ptr == nullptr) {
        ptr = opp_mem::dev_malloc<T>(count, opp_queue);
        opp_consts.push_back((char*)ptr);
    }
}

template <typename T>
T opp_atomic_fetch_add(T* address, T value) {
    return dpct::atomic_fetch_add<sycl::access::address_space::generic_space>(address, value);
}

inline void opp_set_stride(OPP_INT*& data_d, OPP_INT& data_h, OPP_INT new_data) {
    opp_register_const<OPP_INT>(data_d, 1);
    if (data_h != new_data) {
        data_h = new_data;
        opp_mem::copy_host_to_dev<OPP_INT>(data_d, &data_h, 1);
    }    
}