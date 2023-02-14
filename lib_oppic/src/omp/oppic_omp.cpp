
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

//****************************************
void oppic_init(int argc, char **argv, int diags)
{
    oppic_init_core(argc, argv, diags);
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
oppic_dat oppic_decl_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name)
{
    return oppic_decl_dat_core(set, dim, type, size, data, name);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, map_with_cell_index);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, idx, map, acc, map_with_cell_index);
}
oppic_arg oppic_arg_dat(oppic_dat dat, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, acc, map_with_cell_index);
}
oppic_arg oppic_arg_dat(oppic_map map, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(map, acc, map_with_cell_index);
}

//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
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
oppic_dat oppic_decl_particle_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index)
{
    return oppic_decl_particle_dat_core(set, dim, type, size, data, name, cell_index);
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

    oppic_finalize_particle_move_omp(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) printf("oppic_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        oppic_particle_sort(set);
    }
}

//****************************************
void oppic_finalize_particle_move_omp(oppic_set set)
{
    if (OP_DEBUG) printf("oppic_finalize_particle_move_omp set [%s] with particle_remove_count [%d]\n", set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    #pragma omp parallel for
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {
        oppic_dat current_oppic_dat = set->particle_dats->at(i);
        int removed_count = 0;

        for (int j = 0; j < set->size; j++)
        {
            if (set->particle_statuses[j] != NEED_REMOVE) continue;

            char* dat_removed_ptr = (char *)(current_oppic_dat->data + (j * current_oppic_dat->size));
            char* dat_to_replace_ptr = (char *)(current_oppic_dat->data + ((set->size - removed_count - 1) * current_oppic_dat->size));
            
            // Get the last element and replace the hole // Not the Optimum!!!
            // TODO : Can we make NULL data and handle it in sort?
            memcpy(dat_removed_ptr, dat_to_replace_ptr, current_oppic_dat->size); 

            removed_count++;
        }

        // current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, (size_t)(set->size - removed_count) * (size_t)current_oppic_dat->size);
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
    free(set->particle_statuses);
}

//****************************************
void oppic_particle_sort(oppic_set set)
{ TRACE_ME;
    
    if (OP_DEBUG) printf("oppic_particle_sort set [%s]\n", set->name);
    
    int* cell_index_data = (int*)set->cell_index_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(cell_index_data, set->size);

    #pragma omp parallel for
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        auto& dat = set->particle_dats->at(i);
        char *new_data = (char *)malloc(set->array_capacity * dat->size);
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
    oppic_print_dat_to_txtfile_core(dat, file_name_prefix, file_name_suffix);
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