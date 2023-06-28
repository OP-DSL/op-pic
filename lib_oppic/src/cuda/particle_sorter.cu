
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
#include <thrust/generate.h>
#include <thrust/random.h>
#include <thrust/copy.h>

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
            thrust::uniform_int_distribution<int> dist(0, 1000);
            return dist(rng); // Assign random number < 100
        }
    }
};

thrust::device_vector<int> cellIdx_dv;
thrust::device_vector<int> i_dv;

// This assumes all the device data to be valid
void particle_sort_cuda(oppic_set set, bool hole_filling)
{ 
    int set_capacity = set->set_capacity;
    int set_size_plus_removed = set->size + set->particle_remove_count;

    if (OP_DEBUG) 
        opp_printf("particle_sort_cuda", 
        "set [%s] with set capacity [%d] set size+removed [%d] hole_filling [%s]", 
        set->name, set_capacity, set_size_plus_removed, (hole_filling ? "TRUE" : "FALSE"));

    int* cellIdx_dp = (int*)set->mesh_relation_dat->data_d;
    
    cellIdx_dv.reserve(set_size_plus_removed);
    cellIdx_dv.resize(set_size_plus_removed);

    // copy the cell index to the thrust vector
    thrust::copy(cellIdx_dp, cellIdx_dp + set_size_plus_removed, cellIdx_dv.begin());

    if (hole_filling)
    {
        // in hole filling, randomize the cell indices to minimize shared memory issues
        thrust::transform(cellIdx_dv.begin(), cellIdx_dv.end(), 
            cellIdx_dv.begin(), RandomFunctor(seed));
        
        // thrust::replace_if(
        //     cellIdx_dv.begin(), cellIdx_dv.end(),
        //     thrust::placeholders::_1 != MAX_CELL_INDEX,  // Replace all values != MAX_CELL_INDEX
        //     0                                            // New value to assign (zero in this case)
        // );
    }

    i_dv.reserve(set_size_plus_removed);
    i_dv.resize(set_size_plus_removed);
    thrust::sequence(i_dv.begin(), i_dv.end());

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
            sort_dat_according_to_index_int(dat, i_dv, set_capacity, set_size_plus_removed);
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            sort_dat_according_to_index_double(dat, i_dv, set_capacity, set_size_plus_removed);
        }
        else
        {
            std::cerr << "particle_sort_cuda not implemented for data type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());
}

void sort_dat_according_to_index_int(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size)
{ 

    thrust::device_vector<int>* dat_dv = dat->thrust_int;
    thrust::device_vector<int>* sorted_dat_dv = dat->thrust_int_sort;

    copy_according_to_index<int>(dat_dv, sorted_dat_dv, new_idx_dv, set_capacity, size, dat->dim);

    thrust::device_vector<int>* tmp = dat->thrust_int;
    dat->thrust_int = dat->thrust_int_sort;
    dat->thrust_int_sort = tmp;

    dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_int->data());
}

void sort_dat_according_to_index_double(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size)
{ 

    thrust::device_vector<double>* dat_dv = dat->thrust_real;
    thrust::device_vector<double>* sorted_dat_dv = dat->thrust_real_sort;

    copy_according_to_index<double>(dat_dv, sorted_dat_dv, new_idx_dv, set_capacity, size, dat->dim);

    thrust::device_vector<double>* tmp = dat->thrust_real;
    dat->thrust_real = dat->thrust_real_sort;
    dat->thrust_real_sort = tmp;

    dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_real->data());
}