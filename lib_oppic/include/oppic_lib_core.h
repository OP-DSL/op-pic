#pragma once

#include <vector>
#include <map>
#include <iostream>
#include <trace.h>
#include <oppic_util.h>
#include <cstring>

//*************************************************************************************************
#define OP_DEBUG        false

#define OP_READ        0
#define OP_WRITE       1
#define OP_RW          2
#define OP_INC         3
#define OP_MIN         4
#define OP_MAX         5

#define OP_ARG_GBL     0
#define OP_ARG_DAT     1

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

enum MoveStatus 
{
    MOVE_DONE = 1,
    NEED_MOVE,
    NEED_REMOVE,
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
typedef int oppic_arg_type;     /* holds OP_ARG_GBL, OP_ARG_DAT */

struct oppic_arg {
    int index;                  /* index */
    oppic_dat dat;              /* dataset */
    oppic_map map;              /* indirect mapping */
    int dim;                    /* dimension of data */
    int idx;                    /* size (for sequential execution) */
    int size;                   /* size (for sequential execution) */
    char *data;                 /* data on host */
    int *map_data;              /* data on host */
    char *data_d;               /* data on device */
    int *map_data_d;            /* data on device */
    char const *type;           /* datatype */
    oppic_arg_type argtype;
    oppic_access acc;              /* oppic_accessor OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
};

struct oppic_set_core {
    int index;                   /* index */
    int size;                    /* number of elements in set */
    char const *name;            /* name of set */

    bool is_particle = false;
    int diff = 0;                /*number of particles to change*/
    std::vector<int> indexes_to_remove;
    oppic_dat cell_index_dat = nullptr;
    std::vector<oppic_dat> particle_dats;
    std::map<int, part_index> cell_index_v_part_index_map;
    oppic_set cells_set;
};

struct oppic_map_core {
    int index;                  /* index */
    oppic_set from;             /* set pointed from */
    oppic_set to;               /* set pointed to */
    int dim;                    /* dimension of pointer */
    int *map;                   /* array defining pointer */
    char const *name;           /* name of pointer */
    int *map_data_d;            /* device array defining pointer */
};

struct oppic_dat_core {
    int index;                  /* index */
    oppic_set set;              /* set on which data is defined */
    int dim;                    /* dimension of data */
    int size;                   /* size of each element in dataset */
    char *data;                 /* data */
    char const *type;           /* datatype */
    char const *name;           /* name of dataset */
    char *data_d;               /* device data */

    std::vector<char*> thread_data;
    bool is_cell_index = false;
};


extern std::vector<oppic_set> oppic_sets;
extern std::vector<oppic_map> oppic_maps;
extern std::vector<oppic_dat> oppic_dats;


//*************************************************************************************************
// oppic API calls

void oppic_init_core(int argc, char **argv, int diags);
void oppic_exit_core();

oppic_set oppic_decl_set_core(int size, char const *name);

oppic_map oppic_decl_map_core(oppic_set from, oppic_set to, int dim, int *imap, char const *name);

oppic_dat oppic_decl_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name);

oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_dat dat, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat_core(oppic_map map, oppic_access acc, bool map_with_cell_index = false);

// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl_core(double *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl_core(const bool *data, int dim, char const *typ, oppic_access acc);

oppic_set oppic_decl_particle_set_core(int size, char const *name, oppic_set cells_set);
oppic_set oppic_decl_particle_set_core(char const *name, oppic_set cells_set);

oppic_dat oppic_decl_particle_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index = false);

void oppic_increase_particle_count_core(oppic_set particles_set, const int num_particles_to_insert);

void oppic_reset_num_particles_to_insert_core(oppic_set set);

void oppic_mark_particle_to_remove_core(oppic_set set, int particle_index);

void oppic_remove_marked_particles_from_set_core(oppic_set set);
void oppic_remove_marked_particles_from_set_core(oppic_set set, std::vector<int>& idx_to_remove);

void oppic_particle_sort_core(oppic_set set);

//*************************************************************************************************
