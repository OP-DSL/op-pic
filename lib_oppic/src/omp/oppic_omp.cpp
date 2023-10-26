
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

#include <oppic_omp.h>

int OPP_nthreads = 1;
std::vector<int> part_remove_count_per_thr;
std::vector<opp_move_var> move_var_per_thr;

//****************************************
void opp_init(int argc, char **argv)
{
#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULL, "opp::PetscSEQ");
#endif

    oppic_init_core(argc, argv);

    OPP_nthreads = omp_get_max_threads();
}

//****************************************
void opp_exit()
{
#ifdef USE_PETSC
    PetscFinalize();
#endif

    oppic_exit_core();
}

//****************************************
void opp_abort(std::string s)
{
    opp_printf("opp_abort", "%s", s.c_str());
    exit(-1);
}

//****************************************
oppic_set opp_decl_mesh_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map opp_decl_mesh_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    return oppic_decl_map_core(from, to, dim, imap, name);
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
        exit(-1);        
    }
}

//****************************************
bool opp_inc_part_count_with_distribution_omp(oppic_set particles_set, int num_particles_to_insert, oppic_dat part_dist)
{
    opp_profiler->start("IncPartCountWithDistribution");

    if (!oppic_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_profiler->end("IncPartCountWithDistribution");
        return false;
    }

    int* part_mesh_connectivity = (int *)particles_set->mesh_relation_dat->data;
    int* distribution           = (int *)part_dist->data;

    int startX = (particles_set->size - particles_set->diff);
    int distribution_set_size = part_dist->set->size;

    int nthreads = omp_get_max_threads();

    int set_size = particles_set->diff;

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t i = start; i < finish; i++)
        { 
            // iterating for all particles here, makes the performance hit.
            // can we change here?
            for (int j = 0; j < distribution_set_size; j++) 
            {
                if (i < distribution[j]) 
                {
                    part_mesh_connectivity[startX + i] = j;
                    break;
                }  
            }
        }
    } 

    opp_profiler->end("IncPartCountWithDistribution");

    return true;
}

//****************************************
void opp_inc_part_count_with_distribution(oppic_set particles_set, int num_particles_to_insert, oppic_dat part_dist, bool calc_new)
{
    if (OP_DEBUG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    // if (!opp_inc_part_count_with_distribution_omp(particles_set, num_particles_to_insert, part_dist))
    // {
    //     opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_omp failed for particle set [%s]", particles_set->name);
    //     exit(-1);        
    // }

    // perf of core is better than omp :(

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        exit(-1);        
    }
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
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args)
{
    part_remove_count_per_thr.resize(OPP_nthreads);
    std::fill(part_remove_count_per_thr.begin(), part_remove_count_per_thr.end(), 0);

    move_var_per_thr.clear();
    move_var_per_thr.resize(OPP_nthreads);

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    oppic_init_particle_move_core(set);
}

//****************************************
void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
{
    oppic_mark_particle_to_move_core(set, particle_index, move_status);
}

//****************************************
bool opp_finalize_particle_move(oppic_set set)
{ 

    for (int i = 0; i < OPP_nthreads; i++) 
        set->particle_remove_count += part_remove_count_per_thr[i];

    oppic_finalize_particle_move_omp(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) printf("\toppic_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        oppic_particle_sort(set);
    }

    return true;
}

#ifdef USE_OLD_FINALIZE

//****************************************
void oppic_finalize_particle_move_omp(oppic_set set)
{
    if (OP_DEBUG) opp_printf("oppic_finalize_particle_move_omp", "set [%s] with particle_remove_count [%d]\n", 
        set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        // getting a backup of cell index since it will also be rearranged using a random OMP thread
        int *mesh_relation_data = (int *)malloc(set->set_capacity * set->mesh_relation_dat->size); 
        
        memcpy((char*)mesh_relation_data, set->mesh_relation_dat->data, 
            set->set_capacity * set->mesh_relation_dat->size);

        #pragma omp parallel for
        for (int i = 0; i < (int)set->particle_dats->size(); i++)
        {
            oppic_dat current_oppic_dat = set->particle_dats->at(i);
            int removed_count = 0;
            int skip_count = 0;

            for (size_t j = 0; j < set->size; j++)
            {
                if (mesh_relation_data[j] != MAX_CELL_INDEX) continue;

                char* dat_removed_ptr = (char *)(current_oppic_dat->data + (size_t)(j * current_oppic_dat->size));

                // BUG_FIX: (set->size - removed_count - 1) This index could marked to be removed, and if marked, 
                // then there could be an array index out of bounds access error in the future
                while (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX)
                {
                    skip_count++;
                }
                if (j >= (set->size - removed_count - skip_count - 1)) 
                {
                    if (OP_DEBUG) 
                        opp_printf("oppic_finalize_particle_move_omp", 
                            "Current Iteration index [%d] and replacement index %d; hence breaking\n", 
                            j, (set->size - removed_count - skip_count - 1));
                    break;
                }

                size_t offset_byte = (size_t)(set->size - removed_count - skip_count - 1) * current_oppic_dat->size;
                char* dat_to_replace_ptr = (char *)(current_oppic_dat->data + offset_byte);
                
                // Get the last element and replace the hole // Not the Optimum!!!
                // TODO : Can we make NULL data and handle it in sort?
                memcpy(dat_removed_ptr, dat_to_replace_ptr, current_oppic_dat->size); 

                removed_count++;
            }

            // current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, 
            //       (size_t)(set->size - removed_count) * (size_t)current_oppic_dat->size);
        }

        free(mesh_relation_data);
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
}

#else

//****************************************
void oppic_finalize_particle_move_omp(oppic_set set)
{
    if (OP_DEBUG) opp_printf("oppic_finalize_particle_move_omp", "set [%s] with particle_remove_count [%d]\n", 
        set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        int *mesh_relation_data = (int *)set->mesh_relation_dat->data;
        std::vector<std::pair<size_t, size_t>> swap_indices;    // contain hole index and the index from back to swap

        // Idea: The last available element should be copied to the hole
        // In the below scope we try to calculate the element to be swapped with the hole
        {
            // set->particle_remove_count   // the particle count that should be removed
            int removed_count = 0;          // how many elements currently being removed
            int skip_count = 0;             // how many elements from the back is skipped ..
                                            // .. due to that element is also to be removed

            for (size_t j = 0; j < (size_t)set->size; j++)
            {
                // skip if the current index is not to be removed
                if (mesh_relation_data[j] != MAX_CELL_INDEX) 
                    continue;

                // handle if the element from the back is also to be removed
                while ((set->size - removed_count - skip_count - 1 >= 0) && 
                    (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX))
                {
                    skip_count++;
                }

                // check whether the holes are at the back!
                if ((set->size - removed_count - skip_count - 1 < 0) ||
                    (j >= (size_t)(set->size - removed_count - skip_count - 1))) 
                {
                    if (OP_DEBUG) 
                        opp_printf("oppic_finalize_particle_move_core", 
                        "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
                        j, (set->size - removed_count - skip_count - 1), set->name);
                    break;
                }

                swap_indices.push_back(std::make_pair(j, (size_t)(set->size - removed_count - skip_count - 1)));

                removed_count++;
            }
        }

        // For all the dats, fill the holes using the swap_indices
        #pragma omp parallel for
        for (size_t i = 0; i < set->particle_dats->size(); i++)
        {
            oppic_dat dat = set->particle_dats->at(i);

            for (const auto& x : swap_indices)
            {
                char* dat_removed_ptr = (char *)(dat->data + (x.first * dat->size));

                size_t offset_byte = x.second * dat->size;
                char* dat_to_replace_ptr = (char *)(dat->data + offset_byte);
                
                memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 
            }
        }
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
}

#endif

//****************************************
void oppic_particle_sort(oppic_set set)
{ 
    
    if (OP_DEBUG) printf("\toppic_particle_sort set [%s]\n", set->name);
    
    int* mesh_relation_data = (int*)set->mesh_relation_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(mesh_relation_data, set->size);

    #pragma omp parallel for
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        auto& dat = set->particle_dats->at(i);
        char *new_data = (char *)malloc(set->set_capacity * dat->size);
        char *old_data = (char*)dat->data;
        
        for (int j = 0; j < set->size; j++)
        {
            memcpy(new_data + j * dat->size, old_data + idx_before_sort[j] * dat->size, dat->size);
        }

        free(dat->data);
        dat->data = new_data;

        if (dat->is_cell_index)
        { 
            int* cell_index_array = (int*)dat->data;
            int current_cell_index = -1, previous_cell_index = -1;
            std::map<int, part_index>& map = *(set->cell_index_v_part_index_map);
            map.clear();

            for (int j = 0; j < set->size; j++)
            {    
                current_cell_index = cell_index_array[j];
            
                if ((current_cell_index != previous_cell_index) && (current_cell_index >= 0))
                {
                    part_index& pi = map[current_cell_index];
                    pi.start = j;

                    if (previous_cell_index >= 0) map[previous_cell_index].end = (j - 1);
                }
                previous_cell_index = current_cell_index;
            }
            map[previous_cell_index].end = (set->size - 1);
        }
    }
}

//****************************************
void oppic_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_o";
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
void opp_reset_dat(oppic_dat dat, char* val, opp_reset reset)
{
    int set_size = dat->set->size;
    // TODO : reset pragma omp
    for (int i = 0; i < set_size; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
void opp_set_dirtybit(int nargs, oppic_arg *args) 
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
int opp_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args) 
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
void opp_mpi_halo_wait_all(int nargs, oppic_arg *args)
{
    // Nothing to execute here
}

//****************************************
bool opp_part_check_status(opp_move_var& m, int map0idx, oppic_set set, 
    int particle_index, int& remove_count, int thread)
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        part_remove_count_per_thr[thread] += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }

    return true;
}

//****************************************
opp_move_var opp_get_move_var(int thread)
{
    // opp_move_var& move_var = move_var_per_thr[thread];

    // move_var.move_status = OPP_MOVE_DONE;
    // move_var.iteration_one = true;

    opp_move_var move_var;

    return move_var;
}
