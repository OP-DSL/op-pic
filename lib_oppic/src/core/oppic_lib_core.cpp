
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

#include <oppic_lib_core.h>


#ifdef USE_TRACE
    int trace::enabled = 1;
    Trace trace::current = Trace("__TRACE_BASE__");
#endif


//****************************************
std::vector<oppic_set> oppic_sets;
std::vector<oppic_map> oppic_maps;
std::vector<oppic_dat> oppic_dats;

int OP_hybrid_gpu      = 0;
int OP_maps_base_index = 0;
int OP_auto_soa        = 0;
int OP_part_alloc_mult = 1;
int OP_auto_sort       = 1;

opp::Params* opp_params = nullptr;

//****************************************
void oppic_init_core(int argc, char **argv, opp::Params* params) 
{
    oppic_sets.clear();
    oppic_maps.clear();
    oppic_dats.clear();

    opp_params = params;
    
    // these will be overidden by args
    OP_auto_sort = params->get<BOOL>("opp_auto_sort");
    OP_part_alloc_mult = params->get<INT>("opp_allocation_multiple");

    for (int n = 1; n < argc; n++) 
    {
        oppic_set_args_core(argv[n]);
    }
}

//****************************************
void oppic_exit_core() 
{
    for (auto& a : oppic_maps) {
        free(a->map);
        free(a);
    }
    oppic_maps.clear();

    for (auto& a : oppic_dats) {
        free(a->data);
        for (int thr = 0; thr < (int)a->thread_data->size(); thr++) { free(a->thread_data->at(thr)); }
        free(a);
    }
    oppic_dats.clear();

    for (auto& a : oppic_sets) {
        delete a->indexes_to_remove;
        delete a->particle_dats;
        delete a->cell_index_v_part_index_map;
        if (a->particle_statuses) free(a->particle_statuses);
        free(a);
    }
    oppic_sets.clear();

    #ifdef USE_TRACE
        std::string trace_file_name = std::string("trace/trace_") + getTimeStr() + ".csv";
        trace::current.write_profile(trace_file_name);
    #endif
}

//****************************************
void oppic_set_args_core(char *argv) 
{
    char temp[64];
    char *pch;

    pch = strstr(argv, "OPP_ALLOC_MULT=");
    if (pch != NULL) 
    {
        strncpy(temp, pch, 20);
        OP_part_alloc_mult = atoi(temp + 15);
        
        printf("\toppic_set_args_core OP_part_alloc_mult = %d\n", OP_part_alloc_mult);
    }

    pch = strstr(argv, "OPP_AUTO_SORT=");
    if (pch != NULL) 
    {
        strncpy(temp, pch, 20);
        OP_auto_sort = atoi(temp + 14);
        
        printf("\toppic_set_args_core OP_auto_sort = %d\n", OP_auto_sort);
        
        if (!(OP_auto_sort == 1 || OP_auto_sort == 0))
            std::cerr << "OPP_AUTO_SORT should be 0 or 1, Not Auto Sorting" << std::endl;
    }
}

//****************************************
oppic_set oppic_decl_set_core(int size, char const *name) 
{
    oppic_set set          = (oppic_set)malloc(sizeof(oppic_set_core));
    set->index             = oppic_sets.size();
    set->size              = size;
    set->name              = copy_str(name);
    set->core_size         = size;
    set->exec_size         = 0;
    set->nonexec_size      = 0;

    set->is_particle       = false;
    set->set_capacity      = size;
    set->diff              = 0;
    set->mesh_relation_dat = NULL;
    set->cells_set         = NULL;

    set->indexes_to_remove           = new std::vector<int>();
    set->particle_dats               = new std::vector<oppic_dat>();
    set->cell_index_v_part_index_map = new std::map<int, part_index>();
    set->particle_statuses           = NULL;
    set->particle_statuses_d         = NULL;

    oppic_sets.push_back(set);
    return set;
}

//****************************************
oppic_map oppic_decl_map_core(oppic_set from, oppic_set to, int dim, int *imap, char const *name) 
{
    if (from == NULL) 
    {
        printf("\t oppic_decl_map error -- invalid 'from' set for map %s\n", name);
        exit(-1);
    }

    if (to == NULL) 
    {
        printf("\toppic_decl_map error -- invalid 'to' set for map %s\n", name);
        exit(-1);
    }

    if (dim <= 0) 
    {
        printf("\toppic_decl_map error -- negative/zero dimension for map %s\n", name);
        exit(-1);
    }

    oppic_map map = (oppic_map)malloc(sizeof(oppic_map_core));
    map->index    = oppic_maps.size();
    map->from     = from;
    map->to       = to;
    map->dim      = dim;

    map->map      = (int *)malloc((size_t)from->size * (size_t)dim * sizeof(int));
    memcpy(map->map, imap, sizeof(int) * from->size * dim);  

    if (OP_maps_base_index == 1) // convert map to 0 based indexing -- i.e. reduce each map value by 1
    {
        for (int i = 0; i < from->size * dim; i++)
            (map->map[i])--;
    }

    map->map_d        = NULL;
    map->name         = copy_str(name);
    map->user_managed = 1;

    oppic_maps.push_back(map);
    return map;
}

//****************************************
oppic_dat oppic_decl_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name) 
{
    if (set == NULL) {
        printf("\toppic_decl_dat error -- invalid set for data: %s\n", name);
        exit(-1);
    }

    if (dim <= 0) {
        printf("\toppic_decl_dat error -- negative/zero dimension for data: %s\n", name);
        exit(-1);
    }

    oppic_dat dat = (oppic_dat)malloc(sizeof(oppic_dat_core));
    dat->index    = oppic_dats.size();
    dat->set      = set;
    dat->dim      = dim;

    if (set->size > 0)
    {
        size_t bytes = (size_t)dim * (size_t)size * (size_t)(set->size + set->exec_size + set->nonexec_size) * sizeof(char);
        dat->data = (char *)malloc(bytes);
        memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    }
    else
    {
        dat->data = NULL;
    }

    dat->data_d        = NULL;
    dat->name          = copy_str(name);
    dat->type          = copy_str(type);
    dat->size          = dim * size;
    dat->user_managed  = 1;
    dat->mpi_buffer    = NULL;
    dat->buffer_d      = NULL;
    dat->buffer_d_r    = NULL;
    dat->dirty_hd      = Dirty::NotDirty;
    dat->dirtybit      = 1;

    dat->thread_data        = new std::vector<char*>();
    dat->is_cell_index      = false;
    dat->thrust_int         = NULL;
    dat->thrust_real        = NULL;
    dat->thrust_int_sort    = NULL;
    dat->thrust_real_sort   = NULL;

    oppic_dats.push_back(dat);
    return dat;
}

//****************************************
oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping) 
{
    if (dat == nullptr) { std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; oppic_arg arg; return arg; }
    
    return oppic_arg_dat_core(dat, idx, map, dat->dim, dat->type, acc, mapping);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, oppic_access acc, opp_mapping mapping) 
{
    if (dat == nullptr) { std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; oppic_arg arg; return arg; }
    
    return oppic_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, acc, mapping);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping) 
{
    oppic_arg arg;
    arg.index       = -1;
    arg.argtype     = OP_ARG_DAT;

    arg.dat         = dat;
    arg.map         = map;
    arg.dim         = dim;
    arg.idx         = idx;
    
    arg.size        = ((dat == NULL) ? -1 : dat->size);
    arg.data        = ((dat == NULL) ? NULL : dat->data);
    arg.data_d      = ((dat == NULL) ? NULL : dat->data_d);
    arg.map_data    = ((idx == -1 || idx == -2) ? NULL : map->map);
    arg.map_data_d  = ((idx == -1 || idx == -2) ? NULL : map->map_d);

    arg.type        = typ;     // Not used
    arg.acc         = acc;
    arg.opt         = 1;
    arg.sent        = 0;
    
    return arg;
}

oppic_arg oppic_arg_dat_core(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, -1, NULL, acc, mapping);
}

// arg.map has the map, can change to mapping data map if required
oppic_arg oppic_arg_dat_core(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    oppic_arg arg;
    arg.argtype     = OP_ARG_MAP;

    arg.dat         = NULL;
    arg.map         = map;
    arg.dim         = data_map->dim;
    arg.idx         = idx;
    
    arg.size        = data_map->from->size;
    arg.data        = (char*)data_map->map;
    arg.data_d      = (char*)data_map->map_d;
    arg.map_data    = ((map == NULL) ? NULL : map->map);
    arg.map_data_d  = ((map == NULL) ? NULL : map->map_d);

    arg.type        = "int";
    arg.acc         = acc;

    return arg;
}

//****************************************
oppic_arg oppic_arg_gbl_core(double *data, int dim, char const *typ, oppic_access acc)
{
    oppic_arg arg;
    arg.argtype     = OP_ARG_GBL;

    arg.dat         = NULL;
    arg.map         = NULL;
    arg.dim         = dim;
    arg.idx         = -1;
    arg.size        = dim * sizeof(double);
    arg.data        = (char*)data;
    arg.map_data    = NULL;
    arg.type        = typ;
    arg.acc         = acc;
    
    return arg;
}

oppic_arg oppic_arg_gbl_core(int *data, int dim, char const *typ, oppic_access acc)
{
    oppic_arg arg;
    arg.argtype     = OP_ARG_GBL;

    arg.dat         = NULL;
    arg.map         = NULL;
    arg.dim         = dim;
    arg.idx         = -1;
    arg.size        = dim * sizeof(int);
    arg.data        = (char*)data;
    arg.map_data    = NULL;
    arg.type        = typ;
    arg.acc         = acc;
    
    return arg;
}

oppic_arg oppic_arg_gbl_core(const bool *data, int dim, char const *typ, oppic_access acc)
{
    oppic_arg arg;
    arg.argtype     = OP_ARG_GBL;

    arg.dat         = NULL;
    arg.map         = NULL;
    arg.dim         = dim;
    arg.idx         = -1;
    arg.size        = dim * sizeof(bool);
    arg.data        = (char*)data;
    arg.data_d      = NULL;
    arg.map_data_d  = NULL;
    arg.map_data    = NULL;
    arg.type        = typ;
    arg.acc         = acc;
    
    return arg;
}

//****************************************
oppic_set oppic_decl_particle_set_core(int size, char const *name, oppic_set cells_set) // oppic_decl_particle_set should make an oppic_set which changes size during the loop
{   
    oppic_set set = oppic_decl_set_core(size, name);

    set->is_particle = true;
    set->cells_set = cells_set;

    return set;
}

oppic_set oppic_decl_particle_set_core(char const *name, oppic_set cells_set) // oppic_decl_particle_set should make an oppic_set which changes size during the loop
{    
    oppic_set set = oppic_decl_particle_set_core(0, name, cells_set);

    return set;
}

//****************************************
oppic_dat oppic_decl_particle_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index) 
{
    if (set->is_particle != true) 
    {
        printf("\toppic_decl_particle_dat error -- set is not a particle set: %s\n", set->name);
        exit(-1);
    }

    oppic_dat dat = oppic_decl_dat_core(set, dim, type, size, data, name);

    if (cell_index) set->mesh_relation_dat = dat;
    dat->is_cell_index = cell_index;

    set->particle_dats->push_back(dat); 
    return dat;
}

//****************************************
void oppic_increase_particle_count_core(oppic_set particles_set, const int num_particles_to_insert)
{ TRACE_ME;

    if (num_particles_to_insert <= 0) return;

    if (OP_DEBUG) 
        printf("\toppic_increase_particle_count_core set [%s] with size [%d]\n", particles_set->name, num_particles_to_insert);

    int new_particle_set_size = particles_set->size + num_particles_to_insert;

    if (particles_set->set_capacity >= new_particle_set_size)
    {
        if (OP_DEBUG) 
            printf("\toppic_increase_particle_count_core set [%s] No need to reallocate, new size[%d] set_capacity[%d]\n", particles_set->name, new_particle_set_size, particles_set->set_capacity);        
        
        particles_set->size = new_particle_set_size;
        particles_set->diff = num_particles_to_insert;   
        return;
    }

    int new_particle_set_capacity = particles_set->size + num_particles_to_insert * OP_part_alloc_mult;
    // int new_particle_set_capacity = new_particle_set_size;

    for (auto& current_oppic_dat : *(particles_set->particle_dats))
    {
        if (current_oppic_dat->data == NULL) 
        {
            current_oppic_dat->data = (char *)malloc((size_t)(new_particle_set_capacity * current_oppic_dat->size));
        }
        else
        {
            current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, (size_t)(new_particle_set_capacity * current_oppic_dat->size));
        }

        if (current_oppic_dat->is_cell_index)
        {
            int* mesh_rel_array = (int *)current_oppic_dat->data;
            for (int i = particles_set->size; i < new_particle_set_capacity; i++)
                mesh_rel_array[i] = MAX_CELL_INDEX;
        }
        // else
        // {
        //     memset(current_oppic_dat->data + (particles_set->size * current_oppic_dat->size), 0, (new_particle_set_capacity - particles_set->size) * current_oppic_dat->size);
        // }
    }
    
    particles_set->size         = new_particle_set_size;
    particles_set->set_capacity = new_particle_set_capacity;
    particles_set->diff         = num_particles_to_insert;
}

//****************************************
// Not quite necessary - 
// If oppic_reset_num_particles_to_insert isn't called, set->diff will not reset and oppic_par_looppic_inject__ will loop only on the injected particles
// However it will reset upon oppic_increase_particle_count call
void oppic_reset_num_particles_to_insert_core(oppic_set set)
{
    set->diff = 0;
}

//****************************************
void oppic_init_particle_move_core(oppic_set set)
{
    if (OP_DEBUG) printf("\toppic_init_particle_move_core set [%s]\n", set->name);

    // if (set->particle_statuses) free(set->particle_statuses);

    // set->particle_statuses = (int *)malloc(set->size * sizeof(int));
    // memset(set->particle_statuses, 0, set->size * sizeof(int)); // 0 should be MOVE_DONE

    set->particle_remove_count = 0;
}

//****************************************
void oppic_mark_particle_to_move_core(oppic_set set, int particle_index, int move_status)
{
    if (move_status == (int)OPP_NEED_REMOVE) /*outside the mesh*/
    {  
        if (OP_DEBUG) printf("\toppic_mark_particle_to_move_core set [%s] particle_index [%d]\n", set->name, particle_index);

        set->particle_statuses[particle_index] = OPP_NEED_REMOVE;
        (set->particle_remove_count)++;
    }
    else if (move_status != (int)OPP_MOVE_DONE) 
    {
        std::cerr << "oppic_mark_particle_to_move_core Failed to find the cell - Particle Index " << particle_index << std::endl;
    }
}

//****************************************
void oppic_finalize_particle_move_core(oppic_set set)
{
    if (OP_DEBUG) printf("\toppic_finalize_particle_move_core set [%s] size[%d] with particle_remove_count [%d]\n", set->name, set->size, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        int *mesh_relation_data = (int *)malloc(set->set_capacity * set->mesh_relation_dat->size); // getting a backup of cell index since it will also be rearranged using a random OMP thread
        memcpy((char*)mesh_relation_data, set->mesh_relation_dat->data, set->set_capacity * set->mesh_relation_dat->size);

        for (int i = 0; i < (int)set->particle_dats->size(); i++)
        {
            oppic_dat current_oppic_dat = set->particle_dats->at(i);
            int removed_count = 0;
            int skip_count = 0;

            for (int j = 0; j < set->size; j++)
            {
                if (mesh_relation_data[j] != MAX_CELL_INDEX) continue;

                char* dat_removed_ptr = (char *)(current_oppic_dat->data + (j * current_oppic_dat->size));

                // BUG_FIX: (set->size - removed_count - 1) This index could marked to be removed, and if marked, 
                // then there could be an array index out of bounds access error in the future
                while (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX)
                {
                    skip_count++;
                }
                if (j >= (set->size - removed_count - skip_count - 1)) 
                {
                    if (OP_DEBUG) printf("\toppic_finalize_particle_move_core Current Iteration index [%d] and replacement index %d; hence breaking\n", j, (set->size - removed_count - skip_count - 1));
                    break;
                }

                char* dat_to_replace_ptr = (char *)(current_oppic_dat->data + ((set->size - removed_count - skip_count - 1) * current_oppic_dat->size));
                
                // Get the last element and replace the hole // Not the Optimum!!!
                // TODO : Can we make NULL data and handle it in sort?
                memcpy(dat_removed_ptr, dat_to_replace_ptr, current_oppic_dat->size); 

                removed_count++;
            }

            // current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, (size_t)(set->size - removed_count) * (size_t)current_oppic_dat->size);
        }

        free(mesh_relation_data);
    }
    else
    {
        if (OP_DEBUG) printf("\toppic_finalize_particle_move_core Not processing dats since OP_auto_sort = TRUE\n");
    }

    set->size -= set->particle_remove_count;
}

//****************************************
void oppic_mark_particle_to_remove_core(oppic_set set, int particle_index)
{
    set->indexes_to_remove->push_back(particle_index);
}

//****************************************
void oppic_remove_marked_particles_from_set_core(oppic_set set)
{
    oppic_remove_marked_particles_from_set_core(set, *(set->indexes_to_remove));
}

void oppic_remove_marked_particles_from_set_core(oppic_set set, std::vector<int>& idx_to_remove)
{
    int num_particles_to_remove = idx_to_remove.size();

    std::sort(idx_to_remove.begin(), idx_to_remove.end());

    if (OP_DEBUG) printf("\toppic_remove_marked_particles_from_set set [%s] with size [%d]\n", set->name, num_particles_to_remove);

    if (num_particles_to_remove <= 0) return;

    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {
        oppic_dat current_oppic_dat = set->particle_dats->at(i);
        int removed_count = 0;

        for (const int& j : idx_to_remove)
        {
            if (j < 0 || j > set->size) 
            {
                std::cerr << "oppic_remove_marked_particles_from_set set [" << set->name << "] has idx_to_remove [" << j << "]" << std::endl;
                continue;
            }

            char* dat_removed_ptr = (char *)(current_oppic_dat->data + (j * current_oppic_dat->size));
            char* dat_to_replace_ptr = (char *)(current_oppic_dat->data + ((set->size - removed_count - 1) * current_oppic_dat->size));
            
            // Get the last element and replace the hole // Not the Optimum!!!
            // TODO : Can we make NULL data and handle it in sort?
            memcpy(dat_removed_ptr, dat_to_replace_ptr, current_oppic_dat->size); 

            removed_count++;
        }

        current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, (size_t)(set->size - removed_count) * (size_t)current_oppic_dat->size);
    }

    set->size -= num_particles_to_remove;

    idx_to_remove.clear();
}

//****************************************
void oppic_particle_sort_core(oppic_set set)
{ TRACE_ME;
    
    if (OP_DEBUG) printf("\toppic_particle_sort set [%s]\n", set->name);
    
    int* mesh_relation_data = (int*)set->mesh_relation_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(mesh_relation_data, set->set_capacity);

    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        auto& dat = set->particle_dats->at(i);
        char *new_data = (char *)malloc(set->set_capacity * dat->size);
        char *old_data = (char*)dat->data;
        
        for (int j = 0; j < set->set_capacity; j++)
        {
            memcpy(new_data + j * dat->size, old_data + idx_before_sort[j] * dat->size, dat->size);
        }

        free(dat->data);
        dat->data = new_data;
    }

    // this is for double indirections
    // { 
    //     int current_cell_index = -1, previous_cell_index = -1;
    //     std::map<int, part_index>& map = *(set->cell_index_v_part_index_map);
    //     map.clear();

    //     for (int j = 0; j < set->size; j++)
    //     {    
    //         current_cell_index = mesh_relation_data[j];
        
    //         if ((current_cell_index != previous_cell_index) && (current_cell_index >= 0))
    //         {
    //             part_index& pi = map[current_cell_index];
    //             pi.start = j;

    //             if (previous_cell_index >= 0) map[previous_cell_index].end = (j - 1);
    //         }
    //         previous_cell_index = current_cell_index;
    //     }
    //     map[previous_cell_index].end = (set->size - 1);
    // }
}

//****************************************
void oppic_print_dat_to_txtfile_core(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{ TRACE_ME;

    const std::string file_name = std::string("files/") + file_name_prefix + "_" + file_name_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) 
    {
        printf("\toppic_print_dat_to_txtfile_core can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d\n", dat->set->size, dat->dim) < 0) 
    {
        printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
        exit(2);
    }

    for (int i = 0; i < dat->set->size; i++) 
    {
        // fprintf(fp, "%d", i);

        for (int j = 0; j < dat->dim; j++) 
        {
            if (strcmp(dat->type, "double") == 0 ||
                strcmp(dat->type, "double:soa") == 0 ||
                strcmp(dat->type, "double precision") == 0 ||
                strcmp(dat->type, "real(8)") == 0) 
            {
                if (((double *)dat->data)[i * dat->dim + j] == -0.0) 
                { 
                    ((double *)dat->data)[i * dat->dim + j] = +0.0; 
                }

                if (fprintf(fp, " %+2.25lE", ((double *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if (strcmp(dat->type, "float") == 0 ||
                    strcmp(dat->type, "float:soa") == 0 ||
                    strcmp(dat->type, "real(4)") == 0 ||
                    strcmp(dat->type, "real") == 0) 
            {
                if (fprintf(fp, " %+f", ((float *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if (strcmp(dat->type, "int") == 0 ||
                        strcmp(dat->type, "int:soa") == 0 ||
                        strcmp(dat->type, "int(4)") == 0 ||
                        strcmp(dat->type, "integer") == 0 ||
                        strcmp(dat->type, "integer(4)") == 0) 
            {
                if (fprintf(fp, " %+d", ((int *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if ((strcmp(dat->type, "long") == 0) ||
                        (strcmp(dat->type, "long:soa") == 0)) 
            {
                if (fprintf(fp, " %+ld", ((long *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else 
            {
                printf("\toppic_print_dat_to_txtfile_core Unknown type %s, cannot be written to file %s\n", dat->type, file_name.c_str());
                exit(2);
            }
        }

        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

//****************************************
void oppic_dump_dat_core(oppic_dat data) 
{
    fflush(stdout);

    if (data != NULL) 
    {
        for (int i = 0; i < data->set->size; i++) 
        {
            // printf("\t%d", i);

            for (int j = 0; j < data->dim; j++) 
            {
                if (strncmp("double", data->type, 6) == 0)
                {
                    printf("\t %+2.25lE", ((double *)data->data)[i * data->dim + j]);
                } 
                else if (strncmp("real", data->type, 4) == 0) 
                {
                    printf("\t %f", ((float *)data->data)[i * data->dim + j]);
                } 
                else if (strncmp("integer", data->type, 7) == 0) 
                {
                    printf("\t %d", data->data[i * data->dim + j]);
                } 
                else 
                {
                    printf("\toppic_dump_dat_core Unsupported type for dumping %s\n", data->type);
                    exit(0);
                }
            }

            printf("\t\n");
        }
    }

    fflush(stdout);
}

//****************************************
void oppic_print_map_to_txtfile_core(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
{ TRACE_ME;

    const std::string file_name = std::string("files/") + file_name_prefix + "_" + file_name_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) 
    {
        printf("\toppic_print_map_to_txtfile_core can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d\n", map->from->size, map->dim) < 0) 
    {
        printf("\toppic_print_map_to_txtfile_core error writing to %s\n", file_name.c_str());
        exit(2);
    }

    for (int i = 0; i < map->from->size; i++) 
    {
        // fprintf(fp, "%d", i);

        for (int j = 0; j < map->dim; j++) 
        {
            if (fprintf(fp, " %d", ((int *)map->map)[i * map->dim + j]) < 0) 
            {
                printf("\toppic_print_map_to_txtfile_core error writing to %s\n", file_name.c_str());
                exit(2);
            }
        }

        fprintf(fp, "\n");
    }
    
    fclose(fp);
}
 
//****************************************
void* oppic_load_from_file_core(const char* file_name, int set_size, int dim, char const *type, int size)
{
    int fsize = -1, fdim = -1;
    FILE *fp = NULL;
    bool is_error = false;

	if ((fp = fopen(file_name, "r")) == NULL)
	{
		printf("\toppic_load_from_file - Unable to open file %s\n", file_name);
		exit(-1);
	}
	if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2)
	{
		printf("\toppic_load_from_file - error reading file data from %s\n", file_name);
		exit(-1);
	}
    if (fsize < set_size || fdim != dim)
    {
		printf("\toppic_load_from_file - dim and/or set_size issue in file %s\n", file_name);
		exit(-1);        
    }

    void* data = (void *)malloc((size_t)(set_size * dim * size));

    if (strncmp("double", type, 6) == 0)
    {
        double* d_data = (double*)data;

        for (int n = 0; n < set_size; n++)
        {
            switch (dim)
            {
            case 1:
                if (fscanf(fp, " %lf\n", &d_data[n * dim + 0]) != 1) 
                    is_error = true;
                break;
            case 2:
                if (fscanf(fp, " %lf %lf\n", &d_data[n * dim + 0], &d_data[n * dim + 1]) != 2) 
                    is_error = true;
                break;
            case 3:
                if (fscanf(fp, " %lf %lf %lf\n", &d_data[n * dim + 0], &d_data[n * dim + 1], &d_data[n * dim + 2]) != 3) 
                    is_error = true;
                break;
            case 4:
                if (fscanf(fp, " %lf %lf %lf %lf\n", &d_data[n * dim + 0], &d_data[n * dim + 1], &d_data[n * dim + 2], &d_data[n * dim + 3]) != 4) 
                    is_error = true;
                break;    
            case 16:
                if (fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                    &d_data[n * dim + 0], &d_data[n * dim + 1], &d_data[n * dim + 2], &d_data[n * dim + 3], &d_data[n * dim + 4], &d_data[n * dim + 5], &d_data[n * dim + 6], &d_data[n * dim + 7],
                    &d_data[n * dim + 8], &d_data[n * dim + 9], &d_data[n * dim + 10], &d_data[n * dim + 11], &d_data[n * dim + 12], &d_data[n * dim + 13], &d_data[n * dim + 14], &d_data[n * dim + 15]) != 16) 
                    is_error = true;
                break;
            default: is_error = true;
            }
                    
            if (is_error)
            {
                printf("\toppic_load_from_file - error reading from %s at index %d\n", file_name, n);
                free(data);
                exit(-1);
            }
        }
    } 
    else if (strncmp("int", type, 3) == 0) 
    {
        int* i_data = (int*)data;

        for (int n = 0; n < set_size; n++)
        {
            switch (dim)
            {
            case 1:
                if (fscanf(fp, " %d\n", &i_data[n * dim + 0]) != 1) 
                    is_error = true;
                break;
            case 2:
                if (fscanf(fp, " %d %d\n", &i_data[n * dim + 0], &i_data[n * dim + 1]) != 2) 
                    is_error = true;
                break;
            case 3:
                if (fscanf(fp, " %d %d %d\n", &i_data[n * dim + 0], &i_data[n * dim + 1], &i_data[n * dim + 2]) != 3) 
                    is_error = true;
                break;
            case 4:
                if (fscanf(fp, " %d %d %d %d\n", &i_data[n * dim + 0], &i_data[n * dim + 1], &i_data[n * dim + 2], &i_data[n * dim + 3]) != 4) 
                    is_error = true;
                break;    
            default: is_error = true;
            }
                    
            if (is_error)
            {
                printf("\toppic_load_from_file - error reading from %s at index %d\n", file_name, n);
                free(data);
                exit(-1);
            }
        }
    } 
    else 
    {
        printf("\toppic_load_from_file Unsupported type for loading %s\n", type);
        free(data);
        exit(0);
    }

    return data;
}