
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
 
    cellIdx_dv.reserve(sort_size);
    cellIdx_dv.resize(sort_size);

    // copy the cell index to the thrust vector
    int* cellIdx_dp = (int*)set->mesh_relation_dat->data_d;
    thrust::copy(cellIdx_dp + sort_start_index, cellIdx_dp + set_size_plus_removed, cellIdx_dv.begin());

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

    i_dv.reserve(sort_size);
    i_dv.resize(sort_size);
    thrust::sequence(i_dv.begin(), i_dv.end(), sort_start_index);

    thrust::sort_by_key(cellIdx_dv.begin(), cellIdx_dv.end(), i_dv.begin());

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

    cutilSafeCall(cudaDeviceSynchronize());
}

//****************************************
void sort_dat_according_to_index_int(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool hole_filling, int out_start_idx)
{ 

    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat

    thrust::device_vector<int>* dat_dv = dat->thrust_int;
    thrust::device_vector<int>* sorted_dat_dv = dat->thrust_int_sort;

    copy_according_to_index<int>(dat_dv, sorted_dat_dv, new_idx_dv, 
            set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (sorted_dat_dv->begin() + d * set_capacity + out_start_idx), 
                (sorted_dat_dv->begin() + d * set_capacity + out_start_idx + size), 
                (dat_dv->begin() + d * set_capacity + out_start_idx));
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
    
    thrust::device_vector<double>* dat_dv = dat->thrust_real;
    thrust::device_vector<double>* sorted_dat_dv = dat->thrust_real_sort;

    copy_according_to_index<double>(dat_dv, sorted_dat_dv, new_idx_dv, 
            set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (sorted_dat_dv->begin() + d * set_capacity + out_start_idx), 
                (sorted_dat_dv->begin() + d * set_capacity + out_start_idx + size), 
                (dat_dv->begin() + d * set_capacity + out_start_idx));
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