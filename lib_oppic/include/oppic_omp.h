#pragma once

#include <oppic_lib.h>
#include <omp.h>


template <class T> 
void oppic_create_thread_level_data(oppic_arg arg, T init_value)
{
    oppic_dat dat = arg.dat;
    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("oppic_create_thread_level_data template[%d]\n", nthreads);

    if (dat->thread_data.size() <= 0)
    {
        for (int thr = 0; thr < nthreads; thr++)
        {
            char* thr_data = (char *)malloc((size_t)dat->size * (size_t)(dat->set->size) * sizeof(char));;
            dat->thread_data.push_back(thr_data);
        }
    }

    if ((int)dat->thread_data.size() != nthreads)
    {
        std::cerr << "oppic_create_thread_level_data dat [" << dat->name << "] thread_data not properly created [(int)dat->thread_data.size():" << (int)dat->thread_data.size() << " nthreads:" << nthreads << std::endl;
        return;
    }

    for (int thr = 0; thr < nthreads; thr++)
    {
        std::fill_n((T*)(dat->thread_data[thr]), (dat->dim * dat->set->size), init_value);
    }
}

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
                            ((T*)dat->data)[n] += ((T*)dat->thread_data[thr])[n];
                            break;
                        default:
                            std::cerr << "oppic_reduce_thread_level_data dat [" << dat->name << "] acc [" << (int)arg.acc << "] not implemented" << std::endl;
                    }
                }
            }
        }
    }
}



