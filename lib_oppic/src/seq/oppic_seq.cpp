
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

#include <oppic_seq.h>


//****************************************
void oppic_init(int argc, char **argv, opp::Params* params)
{
    oppic_init_core(argc, argv, params);
}

//****************************************
void oppic_exit()
{
    oppic_exit_core();
}

//****************************************
oppic_set oppic_decl_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map oppic_decl_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    return oppic_decl_map_core(from, to, dim, imap, name);
}

//****************************************
oppic_dat oppic_decl_dat(oppic_set set, int dim, opp_data_type dtype, char *data, char const *name)
{
    char* type = nullptr;
    int size = -1;

    if (dtype == OPP_REAL)
    {
        type = "double";
        size = sizeof(double);
    }
    else if (dtype == OPP_INT)
    {
        type = "int";
        size = sizeof(int);       
    }

    return oppic_decl_dat_core(set, dim, type, size, data, name);
}

//****************************************
oppic_map oppic_decl_map_txt(oppic_set from, oppic_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)oppic_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    oppic_map map = oppic_decl_map(from, to, dim, map_data, name);

    free(map_data);

    return map;
}

//****************************************
oppic_dat oppic_decl_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name)
{
    char* type = nullptr;
    int size = -1;

    if (dtype == OPP_REAL)
    {
        type = "double";
        size = sizeof(double);
    }
    else if (dtype == OPP_INT)
    {
        type = "int";
        size = sizeof(int);       
    }

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type, size);

    oppic_dat dat = oppic_decl_dat_core(set, dim, type, size, dat_data, name);

    free(dat_data);

    return dat;
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, mapping);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_dat dat, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, idx, map, acc, mapping);
}


//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg oppic_arg_gbl(int *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg oppic_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
oppic_set oppic_decl_particle_set(int size, char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(size, name, cells_set);
}
oppic_set oppic_decl_particle_set(char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(name, cells_set);
}

//****************************************
oppic_dat oppic_decl_particle_dat(oppic_set set, int dim, opp_data_type dtype, char *data, char const *name, bool cell_index)
{
    char* type = nullptr;
    int size = -1;

    if (dtype == OPP_REAL)
    {
        type = "double";
        size = sizeof(double);
    }
    else if (dtype == OPP_INT)
    {
        type = "int";
        size = sizeof(int);       
    }

    return oppic_decl_particle_dat_core(set, dim, type, size, data, name, cell_index);
}

//****************************************
oppic_dat oppic_decl_particle_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index)
{
    char* type = nullptr;
    int size = -1;

    if (dtype == OPP_REAL)
    {
        type = "double";
        size = sizeof(double);
    }
    else if (dtype == OPP_INT)
    {
        type = "int";
        size = sizeof(int);       
    }

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type, size);

    oppic_dat dat = oppic_decl_particle_dat_core(set, dim, type, size, dat_data, name, cell_index);

    free(dat_data);

    return dat;
}

//****************************************
void oppic_increase_particle_count(oppic_set particles_set, const int num_particles_to_insert)
{
    oppic_increase_particle_count_core(particles_set, num_particles_to_insert);
}

//****************************************
void oppic_reset_num_particles_to_insert(oppic_set set)
{
    oppic_reset_num_particles_to_insert_core(set);
}

//****************************************
void oppic_mark_particle_to_remove(oppic_set set, int particle_index)
{
    oppic_mark_particle_to_remove_core(set, particle_index);
}

//****************************************
void oppic_remove_marked_particles_from_set(oppic_set set)
{
    oppic_remove_marked_particles_from_set_core(set);
}
void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove)
{
    oppic_remove_marked_particles_from_set_core(set, idx_to_remove);
}

//****************************************
void oppic_particle_sort(oppic_set set)
{ TRACE_ME;

    oppic_particle_sort_core(set);
}

//****************************************
void oppic_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_s";
    oppic_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void oppic_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    oppic_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

//****************************************
void oppic_dump_dat(oppic_dat dat)
{
    oppic_dump_dat_core(dat);
}

//****************************************
void oppic_init_particle_move(oppic_set set)
{ TRACE_ME;

    oppic_init_particle_move_core(set);
}

//****************************************
void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
{
    oppic_mark_particle_to_move_core(set, particle_index, move_status);
}

//****************************************
void oppic_finalize_particle_move(oppic_set set)
{ TRACE_ME;

    oppic_finalize_particle_move_core(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) printf("\toppic_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        oppic_particle_sort(set);
    }
}

//****************************************
void oppic_reset_dat(oppic_dat dat, char* val)
{
    int set_size = dat->set->size;

    for (int i = 0; i < set_size; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
void oppic_mpi_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
            (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
            args[n].acc == OP_RW)) 
        {
            args[n].dat->dirty_hd = Dirty::Device;
        }
    }
}

//****************************************
int oppic_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args) 
{
    // for (int n = 0; n < nargs; n++)
    // {
    //     if (args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == Dirty::Host) 
    //     {
    //         op_download_dat(args[n].dat);
    //     }
    // }
    return set->size;
}

//****************************************