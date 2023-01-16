
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


void particle_sort_cuda(oppic_set set)
{ TRACE_ME;

    if (OP_DEBUG) printf("particle_sort_cuda set [%s]\n", set->name);

    int set_size = set->size;

    thrust::device_ptr<int> cellIdx_dp = thrust::device_pointer_cast((int*)set->cell_index_dat->data_d);
    thrust::device_vector<int> cellIdx_dv(cellIdx_dp, cellIdx_dp + set_size);

    thrust::device_vector<int> i_dv(set_size);
    thrust::sequence(i_dv.begin(), i_dv.end());

    thrust::sort_by_key(cellIdx_dv.begin(), cellIdx_dv.end(), i_dv.begin());

    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        auto& dat = set->particle_dats->at(i);

        if (!(strstr(dat->type, ":soa") != NULL || OP_auto_soa || (dat->dim > 1)))
        {
            std::cerr << "particle_sort_cuda not implemented for non SOA data structures [dat " << dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0)
        {
            sort_dat_according_to_index<int>(dat, i_dv, set_size);
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            sort_dat_according_to_index<double>(dat, i_dv, set_size);
        }
        else
        {
            std::cerr << "particle_sort_cuda not implemented for data type " << dat->type << " [dat " << dat->name << "]" << std::endl;
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());
}


// void particle_sort_cuda(oppic_set set)
// { TRACE_ME;

//     if (OP_DEBUG) printf("particle_sort_cuda set [%s]\n", set->name);

//     int set_size = set->size;
// printf("XXXXX set_size %d\n", set_size);
//     thrust::device_ptr<int> cellIdx_dp = thrust::device_pointer_cast((int*)set->cell_index_dat->data_d);
//     thrust::device_vector<int> cellIdx_dv(cellIdx_dp, cellIdx_dp + set_size);

//     thrust::device_vector<int> i_dv(set_size);
//     thrust::sequence(i_dv.begin(), i_dv.end());

// thrust::host_vector<int> h_cellIdx_dv = cellIdx_dv;
// thrust::host_vector<int> h_i_dv = i_dv;

// std::cout << "h_cellIdx_dv: " ;
// thrust::copy_n(h_cellIdx_dv.begin(), set_size, std::ostream_iterator<int>(std::cout, ","));
// std::cout << std::endl << "h_i_dv: " ;
// thrust::copy_n(h_i_dv.begin(), set_size, std::ostream_iterator<int>(std::cout, ","));
// std::cout << std::endl;

//     thrust::sort_by_key(cellIdx_dv.begin(), cellIdx_dv.end(), i_dv.begin());

// thrust::host_vector<int> h_i_dv1 = i_dv;
// std::cout << "i_dv SORTED: " ;
// thrust::copy_n(h_i_dv1.begin(), set_size, std::ostream_iterator<int>(std::cout, ","));
// std::cout << std::endl;

//     for (int i = 0; i < (int)set->particle_dats->size(); i++) // TODO: Remove this -1 to sort cell index
//     {    
//         auto& dat = set->particle_dats->at(i);

//         if (!(strstr(dat->type, ":soa") != NULL || OP_auto_soa || (dat->dim > 1)))
//         {
//             std::cerr << "particle_sort_cuda not implemented for non SOA data structures [dat " << dat->name << "]" << std::endl;
//         }

//         if (strcmp(dat->type, "int") == 0)
//         {
//             sort_dat_according_to_index<int>(dat, i_dv, set_size);
//         }
//         else if (strcmp(dat->type, "double") == 0)
//         {
//             sort_dat_according_to_index<double>(dat, i_dv, set_size);
//         }
//         else
//         {
//             std::cerr << "particle_sort_cuda not implemented for data type " << dat->type << " [dat " << dat->name << "]" << std::endl;
//         }
//     }

// thrust::device_ptr<int> cellIdx_dp1 = thrust::device_pointer_cast((int*)set->cell_index_dat->data_d);
// thrust::device_vector<int> cellIdx_dv1(cellIdx_dp1, cellIdx_dp1 + set_size);
// thrust::host_vector<int> h_i_dv11 = cellIdx_dv1;
// std::cout << "SORTED CELL IDX SORTED: " ;
// thrust::copy_n(h_i_dv11.begin(), set_size, std::ostream_iterator<int>(std::cout, ","));
// std::cout << std::endl;

//     cutilSafeCall(cudaDeviceSynchronize());
// }