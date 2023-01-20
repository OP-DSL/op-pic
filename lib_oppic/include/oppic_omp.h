
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
#include <omp.h>


// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void oppic_create_thread_level_data(oppic_arg arg, T init_value)
{
    oppic_dat dat = arg.dat;
    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("oppic_create_thread_level_data template[%d]\n", nthreads);

    if (dat->thread_data->size() <= 0)
    {
        for (int thr = 0; thr < nthreads; thr++)
        {
            char* thr_data = (char *)malloc((size_t)dat->size * (size_t)(dat->set->size) * sizeof(char));;
            dat->thread_data->push_back(thr_data);
        }
    }

    if ((int)dat->thread_data->size() != nthreads)
    {
        std::cerr << "oppic_create_thread_level_data dat [" << dat->name << "] thread_data not properly created [(int)dat->thread_data.size():" << (int)dat->thread_data->size() << " nthreads:" << nthreads << std::endl;
        return;
    }

    for (int thr = 0; thr < nthreads; thr++)
    {
        std::fill_n((T*)(dat->thread_data->at(thr)), (dat->dim * dat->set->size), init_value);
    }
}

// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void oppic_reduce_thread_level_data(oppic_arg arg)
{
    oppic_dat dat = arg.dat;
    oppic_set set = dat->set;
    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("oppic_reduce_thread_level_data dat [%s] nthreads [%d]\n", dat->name, nthreads);

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = ((dat->dim * set->size)* thr)/nthreads;
            int finish = ((dat->dim * set->size)*(thr+1))/nthreads;
            
            for (int n = start; n < finish; n++)
            {
                for (int thr = 0; thr < nthreads; thr++)
                {
                    switch (arg.acc)
                    {
                        case OP_INC:
                            ((T*)dat->data)[n] += ((T*)dat->thread_data->at(thr))[n];
                            break;
                        default:
                            std::cerr << "oppic_reduce_thread_level_data dat [" << dat->name << "] acc [" << (int)arg.acc << "] not implemented" << std::endl;
                    }
                }
            }
        }
    }
}



