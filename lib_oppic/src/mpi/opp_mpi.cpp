
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

#include <opp_mpi.h>

MPI_Comm OP_MPI_WORLD;
MPI_Comm OP_MPI_GLOBAL;

//*******************************************************************************
void oppic_init(int argc, char **argv, opp::Params* params) 
{
    if (OP_DEBUG) printf("\topp_init\n");

    int flag = 0;
    MPI_Initialized(&flag);
    if (!flag) 
    {
        MPI_Init(&argc, &argv);
    }
    
    OP_MPI_WORLD = MPI_COMM_WORLD;
    OP_MPI_GLOBAL = MPI_COMM_WORLD;
    
    MPI_Comm_rank(OP_MPI_WORLD, &OPP_my_rank);
    MPI_Comm_size(OP_MPI_WORLD, &OPP_comm_size);

    oppic_init_core(argc, argv, params);
}

//*******************************************************************************
void oppic_exit() 
{
    if (OP_DEBUG) printf("\topp_exit\n");

    {   
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        
        // free(set_import_buffer_size);

        // for (int i = 0; i < OP_import_index; i++)
        //     free(OP_import_list[i]);
        // if (OP_import_list)
        //     free(OP_import_list);
        
        // for (int i = 0; i < OP_export_index; i++)
        //     free(OP_export_list[i]);
        // if (OP_export_list)
        //     free(OP_export_list);
    }

    // opp_rt_exit();
    oppic_exit_core();

    int flag = 0;
    MPI_Finalized(&flag);
    if (!flag)
    {
        MPI_Finalize();
    }
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
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_dat_core(set, dim, type.c_str(), size, data, name);
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
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    oppic_dat dat = oppic_decl_dat_core(set, dim, type.c_str(), size, dat_data, name);

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
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_particle_dat_core(set, dim, type.c_str(), size, data, name, cell_index);
}

//****************************************
oppic_dat oppic_decl_particle_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    oppic_dat dat = oppic_decl_particle_dat_core(set, dim, type.c_str(), size, dat_data, name, cell_index);

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

// //****************************************
// void oppic_mark_particle_to_remove(oppic_set set, int particle_index)
// {
//     oppic_mark_particle_to_remove_core(set, particle_index);
// }

// //****************************************
// void oppic_remove_marked_particles_from_set(oppic_set set)
// {
//     oppic_remove_marked_particles_from_set_core(set);
// }
// void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove)
// {
//     oppic_remove_marked_particles_from_set_core(set, idx_to_remove);
// }

//****************************************
void oppic_particle_sort(oppic_set set)
{ TRACE_ME;

    oppic_particle_sort_core(set);
}

// //****************************************
// void oppic_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
// {
//     std::string prefix = std::string(file_name_prefix) + "_s";
//     oppic_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
// }

// //****************************************
// void oppic_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
// {
//     oppic_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
// }

// //****************************************
// void oppic_dump_dat(oppic_dat dat)
// {
//     oppic_dump_dat_core(dat);
// }

//****************************************
void oppic_init_particle_move(oppic_set set)
{ TRACE_ME;

    oppic_init_particle_move_core(set);
}

// //****************************************
// void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
// {
//     oppic_mark_particle_to_move_core(set, particle_index, move_status);
// }

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
    printf("IMPLEMENT HALO EXCHANGE!!!\n");
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

//*******************************************************************************
void opp_partition(op_set prime_set, op_map prime_map, op_dat data)
{
    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(prime_set, prime_map, data);

    opp_desanitize_all_maps();
}

//*******************************************************************************
void opp_partition_core(op_set prime_set, op_map prime_map, op_dat data)
{
    printf("at partition\n");

    if (prime_map != NULL)
    {
        opp_partition_kway(prime_map); // use parmetis kaway partitioning
    }
    else
    {
        std::cerr << "Partitioning prime_map : NULL - UNSUPPORTED Partitioner Specification" << std::endl;
        exit(-1);
    }

    opp_halo_create();

    // set_import_buffer_size = (int *)malloc(OP_set_index * sizeof(int));
    // for (int i = 0; i < OP_set_index; i++)
    //     set_import_buffer_size[i] = 0;

    // sanity check to identify if the partitioning results in ophan elements
    int ctr = 0;
    for (int i = 0; i < prime_map->from->size; i++)
    {
        if (prime_map->map[2 * i] >= prime_map->to->size && prime_map->map[2 * i + 1] >= prime_map->to->size) ctr++;
    }
    printf("Orphan edges: %d\n", ctr);
}

//*******************************************************************************
void opp_sanitize_all_maps()
{
    for (int i = 0; i < oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        for (int n = 0; n < map->from->size; n++)
        {
            for (int d = 1; d < map->from->size; d++)
            {
                if (map->map[n * map->dim + d] < 0)
                {
                    map->map[n * map->dim + d] = map->map[n * map->dim];
                }
            }
        }
    }
}

//*******************************************************************************
void opp_desanitize_all_maps()
{
    for (int i = 0; i < oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        for (int n = 0; n < map->from->size; n++)
        {
            for (int d = 1; d < map->from->size; d++)
            {
                if (map->map[n * map->dim + d] == map->map[n * map->dim])
                {
                    map->map[n * map->dim + d] = -1;;
                }
            }
        }
    }
}
