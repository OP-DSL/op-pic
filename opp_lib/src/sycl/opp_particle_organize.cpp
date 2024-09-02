
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

#include <opp_sycl.h>

constexpr int opp_const_threads_per_block = 192;

dpct::device_vector<int> hf_from_indices_dv;    // holefill from indices - starting in reverse order)
dpct::device_vector<int> hf_sequence_dv;        // sequence of numbers required for hole filling
dpct::device_vector<int> ps_cell_index_dv;      // cell indices required for sorting/swapping
dpct::device_vector<int> ps_swap_indices_dv;    // swap (from mapping) indices required for sorting/swapping

//****************************************
void opp_particle_sort(opp_set set)
{ 
    particle_sort_device(set, false);
}

//****************************************
unsigned int seed = 123; // Seed for random number generator
struct RandomFunctor
{
    mutable dpl::default_engine rng;

    RandomFunctor(unsigned int seed) : rng(seed) {}

    int operator()(int value) const {
        if (value == MAX_CELL_INDEX) {
            return value; // Keep the value of MAX_CELL_INDEX unchanged
        }
        else {
            dpl::uniform_int_distribution<int> dist(0, 1000000000);
            return dist(rng); // Assign random number < 100
        }
    }
};



//****************************************
template <typename T>
void sort_dat_according_to_index(opp_dat dat, const dpct::device_vector<int> &swap_indices, 
    const int set_capacity, const int size, bool shuffle, const int out_start_idx)
{ 
    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat
    
    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<double>(dat->thrust_real, dat->thrust_real_sort, swap_indices, 
    //         set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);

    T *dat_data = (T *)(dat->data_d);
    T *dat_swap_data = (T *)(dat->data_swap_d);
    const int *swap_indices_ptr = (int *)dpct::get_raw_pointer(swap_indices.data());
    const int dat_dim = dat->dim;

    const int nblocks  = (size - 1) / opp_const_threads_per_block + 1;

    // data will be swapped from dat->data_d to dat->data_swap_d according to swap_indices
    opp_queue->submit([&](sycl::handler &cgh) {
        
        cgh.parallel_for(
            sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
            [=](sycl::nd_item<1> item) {
                copy_from<T>(
                    dat_data, dat_swap_data,
                    swap_indices_ptr,
                    set_capacity, set_capacity, 
                    0, out_start_idx,
                    dat_dim, size, 
                    item);
            });
    });

    if (shuffle && OPP_comm_iteration != 0) {
        
        opp_queue->wait();

        // Now the valid data is in dat->data_swap_d, but since it is not the 1st communicator iteration
        // we only swap the newly added indices (since newly added particle count is small compared to 
        // existing particles), simply copy the required indices from dat->data_swap_d to dat->data_d

        for (int d = 0; d < dat->dim; d++)
            std::copy(dpl::execution::make_device_policy(*opp_queue),
                (dat_swap_data + d * set_capacity + out_start_idx),
                (dat_swap_data + d * set_capacity + out_start_idx + size),
                (dat_data + d * set_capacity + out_start_idx));
    }
    else {
        // Now the valid data is in dat->data_swap_d, simply switch the data_swap_d and data_d pointers
        char *tmp = dat->data_d;
        dat->data_swap_d = dat->data_d;
        dat->data_d = tmp;
    }
}

//****************************************
// This assumes all the device data to be valid
void particle_sort_device(opp_set set, bool shuffle)
{ 
    const int set_capacity = set->set_capacity;
    const int set_size_plus_removed = set->size + set->particle_remove_count;
    
    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (shuffle && OPP_comm_iteration != 0) {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    if (OPP_DBG) 
        opp_printf("particle_sort_device", 
            "size[%d] diff[%d] cap[%d] size+rem[%d] shuffle[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (shuffle ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);
    
    opp_profiler->start("PS_Resize");
    ps_cell_index_dv.reserve(set->set_capacity);
    ps_cell_index_dv.resize(sort_size);
    opp_profiler->end("PS_Resize");

    OPP_INT* cellIdx_dp = (OPP_INT*)set->mesh_relation_dat->data_d; // this is cell id data of particle set

    if (shuffle) {
        opp_profiler->start("PS_Shuffle");
        // randomize the cell indices to minimize shared memory issues in later PIC routines
        // The below will create random numbers for each index, and MAX_CELL_INDEX for removed, 
        std::transform(dpl::execution::make_device_policy(*opp_queue),
                    cellIdx_dp + sort_start_index, cellIdx_dp + set_size_plus_removed, 
                    ps_cell_index_dv.begin(),
                    RandomFunctor(seed));
        opp_profiler->end("PS_Shuffle");
    }
    else {
        // copy the cell index to the thrust vector for sorting
        opp_profiler->start("PS_CopyCID");
        std::copy(dpl::execution::make_device_policy(*opp_queue), 
                    cellIdx_dp + sort_start_index, cellIdx_dp + set_size_plus_removed, 
                    ps_cell_index_dv.begin());
        opp_profiler->end("PS_CopyCID");
    }
    
    // Create a sequence of numbers starting from sort_start_index to be used as swap indices
    opp_profiler->start("PS_Sequence");
    ps_swap_indices_dv.reserve(set->set_capacity);
    ps_swap_indices_dv.resize(sort_size);
    dpct::iota(dpl::execution::make_device_policy(*opp_queue),
                ps_swap_indices_dv.begin(), ps_swap_indices_dv.end(), 
                sort_start_index);
    opp_profiler->end("PS_Sequence");

    // sort ps_swap_indices_dv, with the key ps_cell_index_dv. Both keys and values will be sorted
    opp_profiler->start("PS_SortKey");
    dpct::sort(dpl::execution::make_device_policy(*opp_queue),
                ps_cell_index_dv.begin(), ps_cell_index_dv.end(), 
                ps_swap_indices_dv.begin());
    opp_profiler->end("PS_SortKey");

    opp_profiler->start("PS_Dats");
    for (int i = 0; i < (int)set->particle_dats->size(); i++) {    

        opp_dat& dat = set->particle_dats->at(i);

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_sort_device not implemented for non SOA data structures [dat "
                << dat->name << "]" << std::endl;
            opp_abort();
        }

        if (strcmp(dat->type, "int") == 0) {

            sort_dat_according_to_index<OPP_INT>(dat, ps_swap_indices_dv, 
                                    set_capacity, sort_size, shuffle, sort_start_index);
        }
        else if (strcmp(dat->type, "double") == 0) {

            dpct::has_capability_or_fail(opp_queue->get_device(), {sycl::aspect::fp64});
            sort_dat_according_to_index<OPP_REAL>(dat, ps_swap_indices_dv, 
                                    set_capacity, sort_size, shuffle, sort_start_index);
        }
        else {
            std::cerr << "particle_sort_device not implemented for type " << dat->type << " [dat " 
                << dat->name << "]" << std::endl;
            opp_abort();
        }
    }
    opp_profiler->end("PS_Dats");
}

//****************************************
// This assumes all the device data to be valid
void particle_hole_fill_device(opp_set set)
{ 
    const int set_capacity = set->set_capacity;
    const int part_remove_count = set->particle_remove_count;
    const int set_size_plus_removed = set->size + part_remove_count;
    const int nblocks  = (part_remove_count - 1) / opp_const_threads_per_block + 1;

    if (OPP_DBG) 
        opp_printf("particle_hole_fill_device", "remove=%d set_size+removed=%d capacity=%d", 
        part_remove_count, set_size_plus_removed, set_capacity);

    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // sort OPP_remove_particle_indices_d since it can be in shuffled state
    opp_profiler->start("HF_SORT");
    dpl::sort(dpl::execution::make_device_policy(*opp_queue),
        OPP_remove_particle_indices_d,
        OPP_remove_particle_indices_d + part_remove_count);
    opp_profiler->end("HF_SORT");

    // resize hf_sequence_dv and hf_from_indices_dv if required
    if (hf_sequence_dv.capacity() < sort_size) {
        
        hf_sequence_dv.resize(set_capacity);
        dpct::iota(dpl::execution::make_device_policy(*opp_queue),
                   hf_sequence_dv.begin(), hf_sequence_dv.end(), 0);

        hf_from_indices_dv.reserve(set_size_plus_removed);
    }
    hf_from_indices_dv.resize(set_size_plus_removed);

    // get the particle indices in reverse order whose cell index is not MAX_CELL_INDEX
    opp_profiler->start("HF_COPY_IF");
    auto end_iter1 = dpct::copy_if(dpl::execution::make_device_policy(*opp_queue),
        dpl::make_reverse_iterator(hf_sequence_dv.begin() + set_size_plus_removed),
        dpl::make_reverse_iterator(hf_sequence_dv.begin() + sort_start_index),
        dpl::make_reverse_iterator(((OPP_INT*)set->mesh_relation_dat->data_d) + set_size_plus_removed),
        hf_from_indices_dv.begin(), 
        [](int i) { return i != MAX_CELL_INDEX; });
    hf_from_indices_dv.resize(part_remove_count);
    opp_profiler->end("HF_COPY_IF");

    // For all the dats, fill the holes using the swap_indices
    opp_profiler->start("HF_Dats");
    for (opp_dat& dat : *(set->particle_dats)) {

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_hole_fill_device not implemented for non SOA data structures [dat " 
                << dat->name << "]" << std::endl;
            opp_abort();
        }

        if (strcmp(dat->type, "int") == 0) {

            opp_queue->submit([&](sycl::handler &cgh) {
                
                OPP_INT *dat_data = (OPP_INT *)(dat->data_d);
                const int *from_indices = (int *)dpct::get_raw_pointer(hf_from_indices_dv.data());
                const int *remove_indices = OPP_remove_particle_indices_d;
                const int dat_dim = dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
                    [=](sycl::nd_item<1> item) {
                        copy_from_to<OPP_INT>(
                            dat_data, dat_data,
                            from_indices, remove_indices,
                            set_capacity, set_capacity, 
                            0, sort_start_index,
                            dat_dim, part_remove_count, 
                            item);
                    });
            });
        }
        else if (strcmp(dat->type, "double") == 0) {

            dpct::has_capability_or_fail(opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {

                OPP_REAL *dat_data = (OPP_REAL *)(dat->data_d);
                const int *from_indices = (int *)dpct::get_raw_pointer(hf_from_indices_dv.data());
                const int *remove_indices = OPP_remove_particle_indices_d;
                const int dat_dim = dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
                    [=](sycl::nd_item<1> item) { 
                        copy_from_to<OPP_REAL>(
                            dat_data, dat_data,
                            from_indices, remove_indices,
                            set_capacity, set_capacity, 
                            0, sort_start_index,
                            dat_dim, part_remove_count, 
                            item);
                    });
            });
        }
        else {
            std::cerr << "particle_hole_fill_device not implemented for type " << dat->type << " [dat " 
                << dat->name << "]" << std::endl;
            opp_abort();
        }
    }
    opp_queue->wait();
    opp_profiler->end("HF_Dats");
}