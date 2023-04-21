
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
    MPI_Init(&argc, &argv);
    
    if (OP_DEBUG) opp_printf("oppic_init", "");

    OP_MPI_WORLD = MPI_COMM_WORLD;
    OP_MPI_GLOBAL = MPI_COMM_WORLD;
    
    MPI_Comm_rank(OP_MPI_WORLD, &OPP_my_rank);
    MPI_Comm_size(OP_MPI_WORLD, &OPP_comm_size);

    oppic_init_core(argc, argv, params);
}

//*******************************************************************************
void oppic_exit() 
{
    if (OP_DEBUG) opp_printf("oppic_exit", "");

    {   
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
        
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

    oppic_exit_core();

    MPI_Finalize();
}

//****************************************
oppic_set oppic_decl_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map oppic_decl_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    oppic_map map = oppic_decl_map_core(from, to, dim, imap, name);

    opp_printf("oppic_decl_map", OPP_my_rank, " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

    return map;
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

// //****************************************
// void oppic_dump_dat(oppic_dat dat)
// {
//     oppic_dump_dat_core(dat);
// }

//****************************************
void oppic_init_particle_move(oppic_set set)
{ TRACE_ME;

    OPP_comm_iteration = 1;

    oppic_init_particle_move_core(set);
}

// //****************************************
// void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
// {
//     oppic_mark_particle_to_move_core(set, particle_index, move_status);
// }

//****************************************
bool oppic_finalize_particle_move(oppic_set set)
{ TRACE_ME;

    // send the counts and send the particles  
    opp_part_exchange(set);  

    // Can fill the holes here, since the communicated particles will be added at the end
    oppic_finalize_particle_move_core(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) printf("\toppic_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        oppic_particle_sort(set);
    }

    return opp_part_check_all_done(set); 
}

//****************************************
void oppic_reset_dat(oppic_dat dat, char* val, opp_reset reset)
{
    if (!val)
    {
        opp_printf("oppic_reset_dat", "Error: val is NULL");
        return;
    }

    int start = -1;
    int end = -1;

    switch (reset)
    {
        case OPP_Reset_Core:
            start = 0;
            end = dat->set->core_size;
            break;
        case OPP_Reset_Set:
            start = 0;
            end = dat->set->size;
            break;
        case OPP_Reset_ieh:
            start = dat->set->size;
            end = dat->set->size + dat->set->exec_size;
            break;
        case OPP_Reset_inh:
            start = dat->set->size + dat->set->exec_size;
            end = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
            break;
        default:
            opp_printf("oppic_reset_dat", "Error: opp_reset failure");
    }

    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
void oppic_mpi_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        // TODO : Do not include double indirect reductions
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) && 
            (args[n].acc == OP_WRITE || args[n].acc == OP_RW || 
                (args[n].acc == OP_INC && !is_double_indirect_reduction(args[n])))) 
        {
            args[n].dat->dirtybit = 0;
            // args[n].dat->dirty_hd = Dirty::Device;
        }
    }
}

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

    // sanity check to identify if the partitioning results in ophan elements
    int ctr = 0;
    for (int i = 0; i < prime_map->from->size; i++)
    {
        if (prime_map->map[2 * i] >= prime_map->to->size && prime_map->map[2 * i + 1] >= prime_map->to->size) ctr++;
    }
    opp_printf("opp_partition()", "%s Orphans in prime map [%s]: %d", (ctr > 0) ? "Error: " : "", prime_map->name, ctr);

    opp_part_comm_init(); 

    std::vector<std::vector<int>> set_sizes(oppic_sets.size());

    for (oppic_set set : oppic_sets)
    {
        std::vector<int>& recv_vec = set_sizes[set->index];
        recv_vec.resize(OPP_comm_size * 3);

        std::vector<int> sizes{ set->size, set->exec_size, set->nonexec_size };
        MPI_Gather(&(sizes[0]), 3, MPI_INT, &(recv_vec[0]), 3, MPI_INT, OPP_MPI_ROOT, OP_MPI_WORLD);
    }

    // print the set sizes of all ranks after partitioning
    if (OPP_my_rank == OPP_MPI_ROOT)
    {
        std::string log = "";

        for (oppic_set set : oppic_sets)
            log += "\t - " + std::string(set->name);

        opp_printf("opp_partition()", "(size|ieh|inh) %s", log.c_str());

        for (int i = 0; i < OPP_comm_size; i++)
        {
            log = "RANK [" + std::to_string(i) + "]";
            
            for (int j = 0; j < oppic_sets.size(); j++)
                log += "\t\t - " + std::to_string(set_sizes[j][i * 3]) + "|" + 
                    std::to_string(set_sizes[j][i * 3 + 1]) + "|" + std::to_string(set_sizes[j][i * 3 + 2]);

            opp_printf("opp_partition()", "%s", log.c_str());
        }
    }
}

//*******************************************************************************
void opp_sanitize_all_maps()
{
    for (int i = 0; i < oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_sanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        for (int n = 0; n < map->from->size; n++)
        {
            int positive_mapping = -1;
            std::vector<int> index;

            for (int d = 0; d < map->dim; d++)
            {
                if (map->map[n * map->dim + d] < 0)
                {
                    index.push_back(n * map->dim + d);
                }
                else
                {
                    positive_mapping = map->map[n * map->dim + d];
                }
            }

            if (positive_mapping >= 0)
            {
                for (int i : index)
                {
                    map->map[i] = positive_mapping;
                }
            }
            else
            {
                opp_printf("opp_sanitize_all_maps", "Error: No positive mapping found at %d in map: %s", n, map->name);
            }
        }
    }

    // if (OP_DEBUG)
    // {
    //     for (int i = 0; i < oppic_maps.size(); i++)
    //     {
    //         oppic_map map = oppic_maps[i];
            
    //         opp_printf("opp_sanitize_all_maps", OPP_my_rank, " map: %s | from->size: %d | dim: %d", map->name, map->from->size, map->dim);

    //         for (int n = 0; n < map->from->size; n++)
    //         {
    //             for (int d = 1; d < map->dim; d++)
    //             {
    //                 if (map->map[n * map->dim + d] < 0)
    //                 {
    //                     opp_printf("opp_sanitize_all_maps", OPP_my_rank, "Error: map: %s | ptr: %p | negative mapping at index: %d [%d]", map->name, map->map, n, n * map->dim + d);
    //                 }
    //             }
    //         }
    //     }
    // }
}

//*******************************************************************************
void opp_desanitize_all_maps()
{
    for (int i = 0; i < oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        //if (OP_DEBUG) opp_printf("opp_desanitize_all_maps", OPP_my_rank, " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        for (int n = 0; n < map->from->size; n++)
        {
            for (int d = 1; d < map->dim; d++)
            {
                if (map->map[n * map->dim + d] == map->map[n * map->dim])
                {
                    map->map[n * map->dim + d] = -1;
                }
            }
        }
    }
}


















/*******************************************************************************
 * Routine to declare partition information for a given set
 *******************************************************************************/
void decl_partition(op_set set, int *g_index, int *partition) 
{
    part p = (part)malloc(sizeof(part_core));
    p->set = set;
    p->g_index = g_index;
    p->elem_part = partition;
    p->is_partitioned = 0;
    OP_part_list[set->index] = p;
    OP_part_index++;
}

/*******************************************************************************
 * Routine to get partition range on all mpi ranks for all sets
 *******************************************************************************/
void get_part_range(int **part_range, int my_rank, int comm_size, MPI_Comm Comm) 
{
    (void)my_rank;
    for (int s = 0; s < OP_set_index; s++) 
    {
        op_set set = OP_set_list[s];

        int *sizes = (int *)malloc(sizeof(int) * comm_size);
        MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, Comm);

        part_range[set->index] = (int *)malloc(2 * comm_size * sizeof(int));

        int disp = 0;
        for (int i = 0; i < comm_size; i++) 
        {
            part_range[set->index][2 * i] = disp;
            disp = disp + sizes[i] - 1;
            part_range[set->index][2 * i + 1] = disp;
            disp++;
        #ifdef DEBUG
            if (my_rank == OPP_MPI_ROOT && OP_diags > 5)
                printf("range of %10s in rank %d: %d-%d\n", set->name, i,
                    part_range[set->index][2 * i],
                    part_range[set->index][2 * i + 1]);
        #endif
        }
        free(sizes);
    }
}

/*******************************************************************************
 * Routine to get partition (i.e. mpi rank) where global_index is located and
 * its local index
 *******************************************************************************/
int get_partition(int global_index, int *part_range, int *local_index, int comm_size) 
{
    for (int i = 0; i < comm_size; i++) 
    {
        if (global_index >= part_range[2 * i] && global_index <= part_range[2 * i + 1]) 
        {
            *local_index = global_index - part_range[2 * i];
            return i;
        }
    }

    if (OP_DEBUG) 
    {    
        std::string log = "";
        for (int i = 0; i < comm_size; i++) 
        {
            log += std::to_string(i) + "|" + std::to_string(part_range[2 * i]) + "|" + std::to_string(part_range[2 * i + 1]) + "\n";
        }

        opp_printf("get_partition()", OPP_my_rank, "Error: orphan global index %d part_range->\n%s", global_index, log.c_str());
    }

    MPI_Abort(OP_MPI_WORLD, 2);
    return 1;
}

/*******************************************************************************
 * Routine to convert a local index in to a global index
 *******************************************************************************/
int get_global_index(int local_index, int partition, int *part_range, int comm_size) 
{
    (void)comm_size;
    int g_index = part_range[2 * partition] + local_index;
#ifdef DEBUG
    if (g_index > part_range[2 * (comm_size - 1) + 1] && OP_diags > 2)
        printf("Global index larger than set size\n");
#endif
    return g_index;
}

/*******************************************************************************
 * Routine to find the MPI neighbors given a halo list
 *******************************************************************************/
void find_neighbors_set(halo_list List, int *neighbors, int *sizes, int *ranks_size, int my_rank, int comm_size, MPI_Comm Comm) 
{
    int *temp = (int *)malloc(comm_size * sizeof(int));
    int *r_temp = (int *)malloc(comm_size * comm_size * sizeof(int));

    for (int r = 0; r < comm_size * comm_size; r++)
        r_temp[r] = -99;
    for (int r = 0; r < comm_size; r++)
        temp[r] = -99;

    int n = 0;

    for (int r = 0; r < comm_size; r++) 
    {
        if (List->ranks[r] >= 0)
        temp[List->ranks[r]] = List->sizes[r];
    }

    MPI_Allgather(temp, comm_size, MPI_INT, r_temp, comm_size, MPI_INT, Comm);

    for (int i = 0; i < comm_size; i++) 
    {
        if (i != my_rank)
        {
            if (r_temp[i * comm_size + my_rank] > 0) 
            {
                neighbors[n] = i;
                sizes[n] = r_temp[i * comm_size + my_rank];
                n++;
            }
        }
    }
    *ranks_size = n;
    free(temp);
    free(r_temp);
}

bool is_double_indirect_reduction(oppic_arg& arg)
{
    if (arg.argtype == OP_ARG_DAT && arg.mesh_mapping == OPP_Map_from_Mesh_Rel && 
            arg.idx != -1 && arg.acc == OP_INC)
        return true;
    
    return false;
}