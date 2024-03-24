
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

#include <opp_cuda.h>

//****************************************
unsigned int seed = 123; // Seed for random number generator
struct RandomFunctor
{
    mutable thrust::default_random_engine rng;

    RandomFunctor(unsigned int seed) : rng(seed) {}

    __device__ int operator()(int value) const
    {
        if (value == MAX_CELL_INDEX)
        {
            return value; // Keep the value of MAX_CELL_INDEX unchanged
        }
        else
        {
            thrust::uniform_int_distribution<int> dist(0, 1000000000);
            return dist(rng); // Assign random number < 100

            // thrust::normal_distribution<float> dist(50.0f, 10.0f);
            // return static_cast<int>(dist(rng)); // Convert to integer
        }
    }
};

__global__ void copy_int(const int* in_dat_d, int* out_dat_d, const int* indices, int in_stride, 
                                    int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + tid + d * in_stride] = in_dat_d[in_offset + idx + d * out_stride];
        }
    }
}

__global__ void copy_double(const double* in_dat_d, double* out_dat_d, const int* indices, int in_stride, 
                                    int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + tid + d * in_stride] = in_dat_d[in_offset + idx + d * out_stride];
        }
    }
}

//****************************************

thrust::device_vector<int> cellIdx_dv;
thrust::device_vector<int> i_dv;

//****************************************
// This assumes all the device data to be valid
void particle_sort_cuda(oppic_set set, bool hole_filling)
{ 
    int set_capacity = set->set_capacity;
    int set_size_plus_removed = set->size + set->particle_remove_count;
    
    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (hole_filling && OPP_comm_iteration != 0) 
    {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    if (OP_DEBUG) 
        opp_printf("particle_sort_cuda", 
            "setSize[%d] setDiff[%d] setCap[%d] size+rem[%d] fill[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (hole_filling ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);
    
    opp_profiler->start("XResize");
    cellIdx_dv.reserve(set->set_capacity);
    cellIdx_dv.resize(sort_size);
    opp_profiler->end("XResize");

    opp_profiler->start("XCopy");
    // copy the cell index to the thrust vector
    int* cellIdx_dp = (int*)set->mesh_relation_dat->data_d;
    thrust::copy(cellIdx_dp + sort_start_index, cellIdx_dp + set_size_plus_removed, cellIdx_dv.begin());
    opp_profiler->end("XCopy");

    opp_profiler->start("XHoleFill");
    if (hole_filling)
    {
        // in hole filling, randomize the cell indices to minimize shared memory issues
        // The below will create random numbers for each index, and MAX_CELL_INDEX for removed, 
        // ideally this should not be called cell_Idx_dv, better naming would be something like, random ordering
        thrust::transform(cellIdx_dv.begin(), cellIdx_dv.end(), 
            cellIdx_dv.begin(), RandomFunctor(seed));

        // thrust::replace_if(
        //     cellIdx_dv.begin(), cellIdx_dv.end(),
        //     thrust::placeholders::_1 != MAX_CELL_INDEX,  // Replace all values != MAX_CELL_INDEX
        //     0                                            // New value to assign (zero in this case)
        // );
    }
    opp_profiler->end("XHoleFill");

    opp_profiler->start("XSequence");
    i_dv.reserve(set->set_capacity);
    i_dv.resize(sort_size);
    thrust::sequence(i_dv.begin(), i_dv.end(), sort_start_index);
    opp_profiler->end("XSequence");

    // int dis = (int)thrust::distance(i_dv.begin(), i_dv.end());
    // int dis2 = (int)thrust::distance(cellIdx_dv.begin(), cellIdx_dv.end());
    // opp_printf("SORT", "set->size=%d set_size_plus_removed=%d | size %d capacity %d i_dv=%d cellIdx_dv=%d", set->size, set_size_plus_removed, sort_size, set->set_capacity, dis, dis2);

    // opp_profiler->start("XSortKey");
    if (OPP_comm_iteration == 0) opp_profiler->start("XSortKey0");
    else if (OPP_comm_iteration == 1) opp_profiler->start("XSortKey1");
    else opp_profiler->start("XSortKey");
    thrust::sort_by_key(cellIdx_dv.begin(), cellIdx_dv.end(), i_dv.begin());
    // opp_profiler->end("XSortKey");
    if (OPP_comm_iteration == 0) opp_profiler->end("XSortKey0");
    else if (OPP_comm_iteration == 1) opp_profiler->end("XSortKey1");
    else opp_profiler->end("XSortKey");

    opp_profiler->start("XSortDats");
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        oppic_dat& dat = set->particle_dats->at(i);

        if (!(strstr(dat->type, ":soa") != NULL || OP_auto_soa || (dat->dim > 1)))
        {
            std::cerr << "particle_sort_cuda not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0)
        {
            sort_dat_according_to_index_int(dat, i_dv, set_capacity, sort_size, 
                hole_filling, sort_start_index);
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            sort_dat_according_to_index_double(dat, i_dv, set_capacity, sort_size, 
                hole_filling, sort_start_index);
        }
        else
        {
            std::cerr << "particle_sort_cuda not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }
    opp_profiler->end("XSortDats");
}

//****************************************
void sort_dat_according_to_index_int(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool hole_filling, int out_start_idx)
{ 

    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat

    // NOTE: Both commented thrust routine and cuda_kernel function has approx same performance
    copy_according_to_index<int>(dat->thrust_int, dat->thrust_int_sort, new_idx_dv, 
            set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    // const int nblocks  = (size - 1) / 192 + 1;
    // copy_int <<<nblocks, 192>>> (
    //     (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
    //     (int*)thrust::raw_pointer_cast(dat->thrust_int_sort->data()),
    //     (int*)thrust::raw_pointer_cast(new_idx_dv.data()),
    //     set_capacity,
    //     set_capacity,
    //     0,
    //     out_start_idx,
    //     dat->dim,
    //     size);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        cutilSafeCall(cudaDeviceSynchronize());
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_int->begin() + d * set_capacity + out_start_idx));
    }
    else
    {
        thrust::device_vector<int>* tmp = dat->thrust_int;
        dat->thrust_int = dat->thrust_int_sort;
        dat->thrust_int_sort = tmp;

        dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_int->data());
    }
}

//****************************************
void sort_dat_according_to_index_double(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool hole_filling, int out_start_idx)
{ 

    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat
    
    // NOTE: Both commented thrust routine and cuda_kernel function has approx same performance
    copy_according_to_index<double>(dat->thrust_real, dat->thrust_real_sort, new_idx_dv, 
            set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    // const int nblocks  = (size - 1) / 192 + 1;
    // copy_double <<<nblocks, 192>>> (
    //     (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
    //     (double*)thrust::raw_pointer_cast(dat->thrust_real_sort->data()),
    //     (int*)thrust::raw_pointer_cast(new_idx_dv.data()),
    //     set_capacity,
    //     set_capacity,
    //     0,
    //     out_start_idx,
    //     dat->dim,
    //     size);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        cutilSafeCall(cudaDeviceSynchronize());
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_real->begin() + d * set_capacity + out_start_idx));
    }
    else
    {
        thrust::device_vector<double>* tmp = dat->thrust_real;
        dat->thrust_real = dat->thrust_real_sort;
        dat->thrust_real_sort = tmp;

        dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_real->data());
    }
}

//****************************************
//****************************************
__global__ void copy_intFromTo(const int* in_dat_d, int* out_dat_d, const int* from_idx, const int* to_idx, 
                        int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + t_idx + d * in_stride] = in_dat_d[in_offset + f_idx + d * out_stride];
        }
    }
}

__global__ void copy_doubleFromTo(const double* in_dat_d, double* out_dat_d, const int* from_idx, 
            const int* to_idx, int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + t_idx + d * in_stride] = in_dat_d[in_offset + f_idx + d * out_stride];
        }
    }
}
//****************************************


thrust::device_vector<int> to_indices_dv;

//****************************************
// This assumes all the device data to be valid
void particle_hole_fill_cuda(oppic_set set, bool hole_filling)
{ 
    int set_capacity = set->set_capacity;
    int set_size_plus_removed = set->size + set->particle_remove_count;
    
    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (hole_filling && OPP_comm_iteration != 0) 
    {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    // int *mesh_relation_data = (int *)set->mesh_relation_dat->data;
    opp_profiler->start("Z_CID"); // This takes too much time
    thrust::host_vector<int> mesh_relation_data = *(set->mesh_relation_dat->thrust_int);
    opp_profiler->end("Z_CID");

    // if (OP_DEBUG) 
        opp_printf("particle_hole_fill_cuda", 
            "setSize[%d] setDiff[%d] setCap[%d] size+rem[%d] fill[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (hole_filling ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);
    
    to_indices_dv.reserve(set_capacity);
    to_indices_dv.resize(set_size_plus_removed);
    
    opp_profiler->start("Z_SEQ");
    thrust::sequence(to_indices_dv.begin(), to_indices_dv.end());
    opp_profiler->end("Z_SEQ");

    opp_profiler->start("Z_CPYIF");
    auto end_iter = thrust::copy_if(thrust::device, to_indices_dv.begin() + sort_start_index, to_indices_dv.end(), 
        set->mesh_relation_dat->thrust_int->begin(), to_indices_dv.begin(), 
        [] __device__(int i) { return i == MAX_CELL_INDEX; });
    to_indices_dv.resize(thrust::distance(to_indices_dv.begin(), end_iter));
    opp_profiler->end("Z_CPYIF");

    opp_profiler->start("Z_SORT");
    thrust::sort(to_indices_dv.begin(), to_indices_dv.end());
    opp_profiler->end("Z_SORT");

    thrust::host_vector<int> to_indices_hv = to_indices_dv;
    thrust::host_vector<int> from_indices_hv;

    // Idea: The last available element should be copied to the hole
    // In the below scope we try to calculate the element to be swapped with the hole
    {
        int removed_count = 0;          // how many elements currently being removed
        int skip_count = 0;             // how many elements from the back is skipped due to that element is also to be removed

        opp_profiler->start("Z_ARRANGE"); // This takes too much time
        for (size_t j = 0; j < to_indices_hv.size(); j++)
        {
            // handle if the element from the back is also to be removed
            while ((set_size_plus_removed - removed_count - skip_count - 1 >= 0) && 
                (mesh_relation_data[set_size_plus_removed - removed_count - skip_count - 1] == MAX_CELL_INDEX))
            {
                skip_count++;
            }

            // check whether the holes are at the back!
            if ((set_size_plus_removed - removed_count - skip_count - 1 < 0) ||
                (j >= (size_t)(set_size_plus_removed - removed_count - skip_count - 1))) 
            {
                if (OP_DEBUG) 
                    opp_printf("particle_hole_fill_cuda", 
                    "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
                    j, (set_size_plus_removed - removed_count - skip_count - 1), set->name);
                break;
            }

            from_indices_hv.push_back(set_size_plus_removed - removed_count - skip_count - 1);

            removed_count++;
        }
        opp_profiler->end("Z_ARRANGE");
    }

    thrust::device_vector<int> from_indices_dv = from_indices_hv;
    const int nblocks  = ((int)(from_indices_hv.size()) - 1) / 192 + 1;

// opp_printf("particle_hole_fill_cuda", "from_indices_hv size %zu", from_indices_hv.size());

    opp_profiler->start("Z_DATS");
    // For all the dats, fill the holes using the swap_indices
    for (opp_dat& dat : *(set->particle_dats))
    {
        if (!(strstr(dat->type, ":soa") != NULL || OP_auto_soa || (dat->dim > 1)))
        {
            std::cerr << "particle_hole_fill_cuda not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0)
        {
            copy_intFromTo <<<nblocks, 192>>> (
                (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
                (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
                (int*)thrust::raw_pointer_cast(from_indices_dv.data()),
                (int*)thrust::raw_pointer_cast(to_indices_dv.data()),
                set_capacity,
                set_capacity,
                0,
                sort_start_index,
                dat->dim,
                (int)(from_indices_hv.size()));
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            copy_doubleFromTo <<<nblocks, 192>>> (
                (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
                (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
                (int*)thrust::raw_pointer_cast(from_indices_dv.data()),
                (int*)thrust::raw_pointer_cast(to_indices_dv.data()),
                set_capacity,
                set_capacity,
                0,
                sort_start_index,
                dat->dim,
                (int)(from_indices_hv.size()));
        }
        else
        {
            std::cerr << "particle_hole_fill_cuda not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }
    cutilSafeCall(cudaDeviceSynchronize());
    opp_profiler->end("Z_DATS");
}