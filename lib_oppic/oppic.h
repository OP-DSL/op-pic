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

#include <stdlib.h>
#include <list>
#include <map>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include "trace.h"

#define OP_PARTICLES
#define OP_DEBUG  false

//*************************************************************************************************

#define OP_ID NULL

#define OP_READ        0
#define OP_WRITE       1
#define OP_RW          2
#define OP_INC         3
#define OP_MIN         4
#define OP_MAX         5

#define OP_ARG_GBL     0
#define OP_ARG_DAT     1

enum op_iterate_type
{
    OP_ITERATE_ALL = 1,
    OP_ITERATE_INJECTED,
};

enum MoveStatus 
{
    MOVE_DONE = 1,
    NEED_MOVE,
    NEED_REMOVE,
};

//*************************************************************************************************
struct op_set_core;
typedef struct op_set_core *op_set;

struct op_map_core;
typedef struct op_map_core *op_map;

struct op_dat_core;
typedef struct op_dat_core *op_dat;

typedef int op_access;          // holds OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX
typedef int op_arg_type;        // holds OP_ARG_GBL, OP_ARG_DAT

struct part_index {
    int start;
    int end;
};

struct op_arg {
    int index;                  /* index */
    op_dat dat;                 /* dataset */
    op_map map;                 /* indirect mapping */
    int dim;                    /* dimension of data */
    int idx;                    /* size (for sequential execution) */
    int size;                   /* size (for sequential execution) */
    char *data;                 /* data on host */
    int *map_data;              /* data on host */
    char const *type;           /* datatype */
    op_arg_type argtype;
    op_access acc;              /* op_accessor OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
};

struct op_set_core {
    int index;                   /* index */
    int size;                    /* number of elements in set */
    char const *name;            /* name of set */

#ifdef OP_PARTICLES    
    bool is_particle = false;
    int diff = 0;                /*number of particles to change*/
    std::vector<int> indexes_to_remove;
    op_dat cell_index_dat = nullptr;
    std::vector<op_dat> particle_dats;
    std::map<int, part_index> cell_index_v_part_index_map;
    op_set cells_set;
#endif
};

struct op_map_core {
    int index;                  /* index */
    op_set from;                /* set pointed from */
    op_set to;                  /* set pointed to */
    int dim;                    /* dimension of pointer */
    int *map;                   /* array defining pointer */
    char const *name;           /* name of pointer */
};

struct op_dat_core {
    int index;                  /* index */
    op_set set;                 /* set on which data is defined */
    int dim;                    /* dimension of data */
    int size;                   /* size of each element in dataset */
    char *data;                 /* data */
    char const *type;           /* datatype */
    char const *name;           /* name of dataset */

#ifdef OP_PARTICLES
    std::vector<char*> thread_data;
    bool is_cell_index = false;
    std::vector<char*> data_to_insert;    // NOT USED FOR NOW
#endif
};

//*************************************************************************************************
// oppic API calls

void op_init();
void op_exit();
op_set op_decl_set(int size, char const *name);
op_map op_decl_map(op_set from, op_set to, int dim, int *imap, char const *name);
op_dat op_decl_dat(op_set set, int dim, char const *type, int size, char *data, char const *name);
op_arg op_arg_dat(op_dat dat, int idx, op_map map, int dim, const char *typ, op_access acc, bool map_with_cell_index = false) ;
// template <class T> op_arg op_arg_gbl(T *data, int dim, char const *typ, op_access acc);
op_arg op_arg_gbl(double *data, int dim, char const *typ, op_access acc);
op_arg op_arg_gbl(const bool *data, int dim, char const *typ, op_access acc);

#ifdef OP_PARTICLES  
op_set op_decl_particle_set(int size, char const *name, op_set cells_set);
op_set op_decl_particle_set(char const *name, op_set cells_set);
op_dat op_decl_particle_dat(op_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index = false);
void op_increase_particle_count(op_set particles_set, const int num_particles_to_insert);
void op_reset_num_particles_to_insert(op_set set);
void op_insert_to_particle_dat(op_dat dat, const char* data);           // NOT USED FOR NOW
void op_insert_particles_to_set(op_set set);                            // NOT USED FOR NOW
void op_mark_particle_to_remove(op_set set, int particle_index);
void op_remove_marked_particles_from_set(op_set set);
void op_remove_marked_particles_from_set(op_set set, std::vector<int>& idx_to_remove);
void op_particle_sort(op_set set);
op_arg op_arg_dat(op_dat dat, int idx, op_map map, op_access acc, bool map_with_cell_index = false);
op_arg op_arg_dat(op_dat dat, op_access acc, bool map_with_cell_index = false);
op_arg op_arg_dat(op_map map, op_access acc, bool map_with_cell_index = false);

template <class T> void op_create_thread_level_data(op_arg arg, T init_value)
{
    op_dat dat = arg.dat;
    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("op_create_thread_level_data template[%d]\n", nthreads);

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
        std::cerr << "op_create_thread_level_data dat [" << dat->name << "] thread_data not properly created [(int)dat->thread_data.size():" << (int)dat->thread_data.size() << " nthreads:" << nthreads << std::endl;
        return;
    }

    for (int thr = 0; thr < nthreads; thr++)
    {
        std::fill_n((T*)(dat->thread_data[thr]), (dat->dim * dat->set->size), init_value);
    }
}

template <class T> void op_reduce_thread_level_data(op_arg arg)
{
    op_dat dat = arg.dat;
    op_set set = dat->set;
    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("op_reduce_thread_level_data dat [%s] nthreads [%d]\n", dat->name, nthreads);

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
                            std::cerr << "op_reduce_thread_level_data dat [" << dat->name << "] acc [" << (int)arg.acc << "] not implemented" << std::endl;
                    }
                }
            }
        }
    }

}

template <typename... T, typename... OPARG>
void op_par_loop(void (*kernel)(T *...), char const *name, op_set set, op_iterate_type iter_type,
                 OPARG... arguments) {
    printf("op_par_loop %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}

template <typename... T, typename... OPARG>
void op_par_loop_particle(void (*kernel)(T *...), char const *name, op_set set, op_iterate_type iter_type,
                 OPARG... arguments) {
    printf("op_par_loop_particle %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}
 
#endif

std::string getTimeStr();

