
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

#ifdef USE_PETSC
    #include <petscksp.h>
#endif


MPI_Comm OP_MPI_WORLD;
MPI_Comm OP_MPI_GLOBAL;

op_dat op_mpi_get_data(op_dat dat);

//*******************************************************************************
void opp_init(int argc, char **argv) 
{

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULL, "opp::Petsc");
#else
    MPI_Init(&argc, &argv);
#endif

    OP_MPI_WORLD = MPI_COMM_WORLD;
    OP_MPI_GLOBAL = MPI_COMM_WORLD;
    
    MPI_Comm_rank(OP_MPI_WORLD, &OPP_rank);
    MPI_Comm_size(OP_MPI_WORLD, &OPP_comm_size);

    if (OP_DEBUG) opp_printf("oppic_init", "");
    
    oppic_init_core(argc, argv);
}

//*******************************************************************************
void opp_exit() 
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

#ifdef USE_PETSC
    PetscFinalize();
#else
    MPI_Finalize();
#endif
    
}

//****************************************
void opp_abort()
{
    MPI_Abort(OP_MPI_WORLD, 2);
}

//****************************************
oppic_set opp_decl_mesh_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map opp_decl_mesh_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    oppic_map map = oppic_decl_map_core(from, to, dim, imap, name);

    if (OP_DEBUG) opp_printf("oppic_decl_map", OPP_rank, " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

    return map;
}

//****************************************
oppic_dat opp_decl_mesh_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);
}

//****************************************
oppic_map oppic_decl_map_txt(oppic_set from, oppic_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)oppic_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    oppic_map map = opp_decl_mesh_map(from, to, dim, map_data, name);

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
oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, mapping);
}

//****************************************
oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, acc, mapping);
}
oppic_arg opp_get_arg(oppic_dat dat, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, acc, mapping);
}
oppic_arg opp_get_arg(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, acc, mapping);
}
oppic_arg opp_get_arg(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, idx, map, acc, mapping);
}


//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg opp_get_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg opp_get_arg_gbl(int *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg opp_get_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
oppic_set opp_decl_part_set(int size, char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(size, name, cells_set);
}
oppic_set opp_decl_part_set(char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(name, cells_set);
}

//****************************************
oppic_dat opp_decl_part_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_particle_dat_core(set, dim, type.c_str(), size, (char*)data, name, cell_index);
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
    if (!oppic_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_printf("oppic_increase_particle_count", "Error: oppic_increase_particle_count_core failed for particle set [%s]", particles_set->name);
        MPI_Abort(OP_MPI_WORLD, 1);        
    }
}

//****************************************
void opp_inc_part_count_with_distribution(oppic_set particles_set, int num_particles_to_insert, oppic_dat part_dist)
{
    if (OP_DEBUG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        MPI_Abort(OP_MPI_WORLD, 1);        
    }
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
{ 

    oppic_particle_sort_core(set);
}

//****************************************
void opp_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_m" + std::to_string(OPP_comm_size);
    oppic_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    oppic_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

// //****************************************
// void oppic_dump_dat(oppic_dat dat)
// {
//     oppic_dump_dat_core(dat);
// }

//****************************************
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args)
{ 

    oppic_init_particle_move_core(set);

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }
    else
    {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++)
        {
            if (args[i].argtype == OP_ARG_DAT && args[i].dat->set->is_particle)
            {
                args[i].data = args[i].dat->data;
            }
        }
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
}

// //****************************************
// void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
// {
//     oppic_mark_particle_to_move_core(set, particle_index, move_status);
// }

//****************************************
bool opp_finalize_particle_move(oppic_set set)
{ 

    if (OP_DEBUG) opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    // send the counts and send the particles  
    opp_part_exchange(set);  

    // Can fill the holes here, since the communicated particles will be added at the end
    oppic_finalize_particle_move_core(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) opp_printf("opp_finalize_particle_move", "auto sorting particle set [%s]", set->name);
        oppic_particle_sort(set);
    }

    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop

        return false; // all mpi ranks do not have anything to communicate to any rank
    }
        
    opp_part_wait_all(set); // wait till all the particles are communicated and added to the dats

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    return true;
}

//****************************************
void opp_reset_dat(oppic_dat dat, char* val, opp_reset reset)
{
    if (!val)
    {
        opp_printf("opp_reset_dat", "Error: val is NULL");
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
        case OPP_Reset_All:
            start = 0;
            end = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
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
            opp_printf("opp_reset_dat", "Error: opp_reset failure");
    }

    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
void opp_mpi_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        // TODO : Do not include double indirect reductions
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) && 
            (args[n].acc == OP_WRITE || args[n].acc == OP_RW || 
                (args[n].acc == OP_INC)))  //  && !is_double_indirect_reduction(args[n])
        {
            args[n].dat->dirtybit = 1;
            // args[n].dat->dirty_hd = Dirty::Device;
        }
    }
}

//*******************************************************************************
void opp_partition(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();
}

//*******************************************************************************
void opp_partition_core(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    if (lib_name == "PARMETIS_KWAY")
    {
        if (prime_map != NULL)
        {
            opp_partition_kway(prime_map); // use parmetis kway partitioning
        }
        else
        {
            opp_printf("opp_partition", "Error: Partitioning prime_map : NULL - UNSUPPORTED Partitioner Specification");
            MPI_Abort(OP_MPI_WORLD, 1);
        }
    }
    else if (lib_name == "PARMETIS_GEOM")
    {
        if (data != NULL)
        {
            opp_partition_geom(data); // use parmetis geometric partitioning
        }
        else
        {
            opp_printf("opp_partition", "Error: Partitioning geom dat : NULL - UNSUPPORTED Partitioner Specification");
            MPI_Abort(OP_MPI_WORLD, 1);
        }
    }
    else if (lib_name == "EXTERNAL")
    {
        if (data != NULL)
        {
            opp_partition_external(prime_set, data); // use external partitioning dat
        }
        else
        {
            opp_printf("opp_partition", "Error: Partitioning color dat : NULL - UNSUPPORTED Partitioner Specification");
            MPI_Abort(OP_MPI_WORLD, 1);
        }
    }
    else if (lib_name != "")
    {
        opp_printf("opp_partition", "Error: Unsupported lib_name [%s] - UNSUPPORTED Partitioner Specification", lib_name.c_str());
        MPI_Abort(OP_MPI_WORLD, 1);
    }

    opp_halo_create();

    // TODO : I think this sanity check is wrong, in cell->nodes mapping there can be nodes indexed in import non-exec halo,
    // hence the mapping is obvously greater than to_set->size 
        // sanity check to identify if the partitioning results in ophan elements
        // int ctr = 0;
        // for (int i = 0; i < prime_map->from->size; i++)
        // {
        //     if (prime_map->map[2 * i] >= prime_map->to->size && prime_map->map[2 * i + 1] >= prime_map->to->size) 
        //         ctr++;
        // }
        // opp_printf("opp_partition()", "%s Orphans in prime map [%s]: %d", (ctr > 0) ? "Error:" : "", prime_map->name, ctr);

    opp_part_comm_init(); 

    std::vector<std::vector<int>> set_sizes(oppic_sets.size());

    for (oppic_set set : oppic_sets)
    {
        std::vector<int>& recv_vec = set_sizes[set->index];
        recv_vec.resize(OPP_comm_size * 3);

        std::vector<int> sizes{ set->size, set->exec_size, set->nonexec_size };
        MPI_Gather(&(sizes[0]), 3, MPI_INT, &(recv_vec[0]), 3, MPI_INT, OPP_ROOT, OP_MPI_WORLD);
    }

    // print the set sizes of all ranks after partitioning
    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "";

        for (oppic_set set : oppic_sets)
            log += "\t - " + std::string(set->name);

        opp_printf("opp_partition()", "(size|ieh|inh) %s", log.c_str());

        for (int i = 0; i < OPP_comm_size; i++)
        {
            log = "RANK [" + std::to_string(i) + "]";
            
            for (int j = 0; j < oppic_sets.size(); j++)
                log += "\t- " + std::to_string(set_sizes[j][i * 3]) + "|" + 
                    std::to_string(set_sizes[j][i * 3 + 1]) + "|" + std::to_string(set_sizes[j][i * 3 + 2]);

            opp_printf("opp_partition()", "%s", log.c_str());
        }
    }
}

std::map<int, oppic_dat> negative_mapping_indices;

//*******************************************************************************
void opp_sanitize_all_maps()
{
    for (int i = 0; i < oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_sanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        std::string name = std::string("AUTO_DAT_") + map->name;
        oppic_dat dat = opp_decl_mesh_dat(map->from, map->dim, DT_INT, (char*)map->map, name.c_str());  
        negative_mapping_indices[map->index] = dat;

        memset(dat->data, 0, (map->from->size * dat->size));

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
                    ((int*)dat->data)[i] = 1;
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
            
    //         opp_printf("opp_sanitize_all_maps", OPP_rank, " map: %s | from->size: %d | dim: %d", map->name, map->from->size, map->dim);

    //         for (int n = 0; n < map->from->size; n++)
    //         {
    //             for (int d = 1; d < map->dim; d++)
    //             {
    //                 if (map->map[n * map->dim + d] < 0)
    //                 {
    //                     opp_printf("opp_sanitize_all_maps", OPP_rank, "Error: map: %s | ptr: %p | negative mapping at index: %d [%d]", map->name, map->map, n, n * map->dim + d);
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

        if (OP_DEBUG) opp_printf("opp_desanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        auto it = negative_mapping_indices.find(map->index);
        if (it == negative_mapping_indices.end())
        {
            opp_printf("opp_desanitize_all_maps", "Error: Negative mappings not found for map: %s", map->name);
            continue;
        }
            
        oppic_dat dat = it->second;

        for (int x = 0; x < (map->from->size + map->from->exec_size) * map->dim; x++)
        {
            if (((int*)dat->data)[x] == 1)
                map->map[x] = -1;
        }
    }

    // could realloc the dat to a lower size, if required
    negative_mapping_indices.clear();
}


void opp_mpi_print_dat_to_txtfile(op_dat dat, const char *file_name) 
{
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + std::to_string(OPP_comm_size) + std::string("_") + file_name;
    // rearrange data back to original order in mpi
    op_dat temp = op_mpi_get_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    free(temp->data);
    free(temp->set);
    free(temp);
}

opp_move_var opp_get_move_var()
{
// TODO_IMM : use a buffered opp_move_var instead, could use a global variable and reset

    opp_move_var m;

    if (OPP_comm_iteration != 0) // TRUE means communicated particles, no need to do the iteration one calculations
        m.OPP_iteration_one = false;
    
    return m;
}