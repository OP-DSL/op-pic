
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

#include <vector>
#include <map>
#include <iostream>
#include <trace.h>
#include <opp_params.h>
#include <oppic_util.h>
#include <cstring>
#include <limits.h>

//*************************************************************************************************
#define OP_DEBUG       false

#define MAX_CELL_INDEX     INT_MAX

#define OP_READ        0
#define OP_WRITE       1
#define OP_RW          2
#define OP_INC         3
#define OP_MIN         4
#define OP_MAX         5

#define OP_ARG_GBL     0
#define OP_ARG_DAT     1
#define OP_ARG_MAP     2

#define ZERO_double    0.0
#define ZERO_float     0.0f
#define ZERO_int       0
#define ZERO_uint      0
#define ZERO_ll        0
#define ZERO_ull       0
#define ZERO_bool      0

//*************************************************************************************************
enum oppic_iterate_type
{
    OP_ITERATE_ALL = 1,
    OP_ITERATE_INJECTED,
};

enum move_status 
{
    OPP_MOVE_DONE = 0,
    OPP_NEED_MOVE,
    OPP_NEED_REMOVE,
};

enum DeviceType
{
    Device_CPU = 1,
    Device_GPU = 2,
};

enum Dirty
{
    NotDirty = 0,
    Device = 1,
    Host = 2,
};

struct part_index {
    int start;
    int end;
};

//*************************************************************************************************
struct oppic_set_core;
typedef struct oppic_set_core *oppic_set;

struct oppic_map_core;
typedef struct oppic_map_core *oppic_map;

struct oppic_dat_core;
typedef struct oppic_dat_core *oppic_dat;

typedef int oppic_access;       /* holds OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
typedef int oppic_arg_type;     /* holds OP_ARG_GBL, OP_ARG_DAT, OP_ARG_MAP */

struct oppic_arg {
    int index;                  /* index */
    oppic_dat dat;              /* dataset */
    oppic_map map;              /* indirect mapping */
    int dim;                    /* dimension of data */
    int idx;                    /* size (for sequential execution) */
    int size;                   /* size (for sequential execution) */
    char *data;                 /* data on host */
    char *data_d;               /* data on device (for CUDA execution) */
    int *map_data;              /* data on host */
    int *map_data_d;            /* data on device (for CUDA execution) */
    char const *type;           /* datatype */
    oppic_access acc;           /* oppic_accessor OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
    oppic_arg_type argtype;
    int sent;                   /* flag to indicate if this argument has data in flight under non-blocking MPI comms*/
    int opt;                    /* flag to indicate if this argument is in use */
};

struct oppic_set_core {
    int index;                              /* index */
    int size;                               /* number of elements in set */
    char const *name;                       /* name of set */
    int core_size;                          /* number of core elements in an mpi process*/
    int exec_size;                          /* number of additional imported elements to be executed */
    int nonexec_size;                       /* number of additional imported elements that are not executed */

    bool is_particle;                       /* is it a particle set */
    int array_capacity;                     /* capacity of the allocated array */
    int diff;                               /* number of particles to change */
    std::vector<int>* indexes_to_remove;
    oppic_dat cell_index_dat = NULL;
    std::vector<oppic_dat>* particle_dats;
    std::map<int, part_index>* cell_index_v_part_index_map;
    int* particle_statuses;
    int* particle_statuses_d;
    int particle_remove_count;
    oppic_set cells_set;
};

struct oppic_map_core {
    int index;                  /* index */
    oppic_set from;             /* set pointed from */
    oppic_set to;               /* set pointed to */
    int dim;                    /* dimension of pointer */
    int *map;                   /* array defining pointer */
    int *map_d;                 /* device array defining pointer */
    char const *name;           /* name of pointer */
    int user_managed;           /* indicates whether the user is managing memory */
};

struct oppic_dat_core {
    int index;                  /* index */
    oppic_set set;              /* set on which data is defined */
    int dim;                    /* dimension of data */
    int size;                   /* size of each element in dataset */
    char *data;                 /* data on host */
    char *data_d;               /* data on device (GPU) */
    char const *type;           /* datatype */
    char const *name;           /* name of dataset */
    char *buffer_d;             /* buffer for MPI halo sends on the devidce */
    char *buffer_d_r;           /* buffer for MPI halo receives on the devidce */
    int dirtybit;               /* flag to indicate MPI halo exchange is needed*/
    Dirty dirty_hd;             /* flag to indicate dirty status on host and device */
    int user_managed;           /* indicates whether the user is managing memory */
    void *mpi_buffer;           /* ponter to hold the mpi buffer struct for the op_dat*/    

    std::vector<char*>* thread_data;
    bool is_cell_index;
};

//*************************************************************************************************
// oppic API calls

void oppic_init_core(int argc, char **argv, opp::Params* params);
void oppic_exit_core();

void oppic_set_args_core(char *argv);

oppic_set oppic_decl_set_core(int size, char const *name);

oppic_map oppic_decl_map_core(oppic_set from, oppic_set to, int dim, int *imap, char const *name);

oppic_dat oppic_decl_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name);

oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_dat dat, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_map data_map, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_map data_map, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index = false);

// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl_core(double *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl_core(int *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl_core(const bool *data, int dim, char const *typ, oppic_access acc);

oppic_set oppic_decl_particle_set_core(int size, char const *name, oppic_set cells_set);
oppic_set oppic_decl_particle_set_core(char const *name, oppic_set cells_set);

oppic_dat oppic_decl_particle_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index = false);

void oppic_increase_particle_count_core(oppic_set particles_set, const int num_particles_to_insert);

void oppic_reset_num_particles_to_insert_core(oppic_set set);

void oppic_init_particle_move_core(oppic_set set);

void oppic_mark_particle_to_move_core(oppic_set set, int particle_index, int move_status);

void oppic_finalize_particle_move_core(oppic_set set);

void oppic_mark_particle_to_remove_core(oppic_set set, int particle_index);

void oppic_remove_marked_particles_from_set_core(oppic_set set);
void oppic_remove_marked_particles_from_set_core(oppic_set set, std::vector<int>& idx_to_remove);

void oppic_particle_sort_core(oppic_set set);

void oppic_print_dat_to_txtfile_core(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix);
void oppic_print_map_to_txtfile_core(oppic_map map, const char *file_name_prefix, const char *file_name_suffix);

void oppic_dump_dat_core(oppic_dat data);
//*************************************************************************************************

extern int OP_hybrid_gpu;
extern int OP_maps_base_index;
extern int OP_auto_soa;
extern int OP_part_alloc_mult;
extern int OP_auto_sort;

extern std::vector<oppic_set> oppic_sets;
extern std::vector<oppic_map> oppic_maps;
extern std::vector<oppic_dat> oppic_dats;

extern opp::Params* opp_params;

void* oppic_load_from_file_core(const char* file_name, int set_size, int dim, char const *type, int size);