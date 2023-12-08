
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

#ifdef USE_MPI
    #include <opp_mpi.h>
#endif

using namespace opp;

//****************************************
std::vector<oppic_set> oppic_sets;
std::vector<oppic_map> oppic_maps;
std::vector<oppic_dat> oppic_dats;

int OP_hybrid_gpu                   = 0;
int OP_maps_base_index              = 0;
int OP_auto_soa                     = 0;
int OP_part_alloc_mult              = 1;
int OP_auto_sort                    = 1;
int OPP_mpi_part_alloc_mult         = 1;
int OPP_rank                        = 0;
int OPP_comm_size                   = 1;
int OPP_comm_iteration              = 0;
int OPP_max_comm_iteration          = 0;
int OPP_iter_start                  = 0;
int OPP_iter_end                    = 0;
int OPP_main_loop_iter              = 0;
int OPP_gpu_threads_per_block       = OPP_DEFAULT_GPU_THREADS_PER_BLOCK;
size_t OPP_gpu_shared_mem_per_block = -1;
int *OPP_mesh_relation_data         = nullptr;
int *OPP_mesh_relation_data_d       = nullptr;
int OPP_part_cells_set_size         = 0;

std::unique_ptr<opp::Params> opp_params;
std::unique_ptr<opp::Profiler> opp_profiler;

std::shared_ptr<BoundingBox> boundingBox;
std::shared_ptr<CellMapper> cellMapper;
std::shared_ptr<Comm> comm;
std::unique_ptr<GlobalParticleMover> globalMover;
bool useGlobalMove = true;

//****************************************
void oppic_init_core(int argc, char **argv) 
{
    oppic_sets.clear();
    oppic_maps.clear();
    oppic_dats.clear();

    opp_params = std::make_unique<opp::Params>(argv[1]);
    opp_profiler = std::make_unique<opp::Profiler>();

    // these will be overidden by args
    OP_auto_sort = opp_params->get<OPP_BOOL>("opp_auto_sort");
    OP_part_alloc_mult = opp_params->get<OPP_INT>("opp_allocation_multiple");

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
        free((char*)a->name);
        free(a);
    }
    oppic_maps.clear();

    for (auto& a : oppic_dats) {
        free(a->data);
        for (int thr = 1; thr < (int)a->thread_data->size(); thr++) { 
            free(a->thread_data->at(thr)); 
        }
        delete a->thread_data;
        free((char*)a->name);
        free((char*)a->type);
        free(a);
    }
    oppic_dats.clear();

    for (auto& a : oppic_sets) {
        delete a->indexes_to_remove;
        delete a->particle_dats;
        delete a->cell_index_v_part_index_map;
        if (a->particle_statuses) free(a->particle_statuses);
        free((char*)a->name);
        free(a);
    }
    oppic_sets.clear();

    if (opp_profiler.get())
    {    
        if (opp_params->get<OPP_BOOL>("opp_profile_all"))
            opp_profiler->printProfile(true);
        opp_profiler->printProfile();
    }
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

    pch = strstr(argv, "OPP_MPI_ALLOC_MULT=");
    if (pch != NULL) 
    {
        strncpy(temp, pch, 20);
        OPP_mpi_part_alloc_mult = atoi(temp + 15);
        
        printf("\toppic_set_args_core OPP_mpi_part_alloc_mult = %d\n", OPP_mpi_part_alloc_mult);
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
    set->particle_size     = 0;

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
        opp_printf("oppic_decl_map", "Error -- invalid 'from' set for map %s", name);
        exit(-1);
    }

    if (to == NULL) 
    {
        opp_printf("oppic_decl_map", "Error -- invalid 'to' set for map %s", name);
        exit(-1);
    }

    if (dim <= 0) 
    {
        opp_printf("oppic_decl_map", "Error -- negative/zero dimension for map %s", name);
        exit(-1);
    }

    if (from->is_particle || to->is_particle) 
    {
        opp_printf("oppic_decl_map", 
            "Error -- cannot have mappings between a dynamic (particle) set [%s to %s]", 
            from->name, to->name);
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
    if (set == NULL) 
    {
        printf("\toppic_decl_dat error -- invalid set for data: %s\n", name);
        exit(-1);
    }

    if (dim <= 0) 
    {
        printf("\toppic_decl_dat error -- negative/zero dimension for data: %s\n", name);
        exit(-1);
    }

    oppic_dat dat = (oppic_dat)malloc(sizeof(oppic_dat_core));
    dat->index    = oppic_dats.size();
    dat->set      = set;
    dat->dim      = dim;

    if (set->size > 0)
    {
        size_t bytes = (size_t)dim * (size_t)size * 
                            (size_t)(set->size + set->exec_size + set->nonexec_size) * sizeof(char);
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
    dat->mpi_reduc_buffer = NULL;
    dat->reduc_comm    = OPP_Reduc_NO_Comm;
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

    set->particle_size += dat->size;

    oppic_dats.push_back(dat);
    return dat;
}

//****************************************
oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping) 
{
    if (dat == nullptr) 
    { 
        std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; 
        oppic_arg arg; 
        return arg; 
    }
    
    return oppic_arg_dat_core(dat, idx, map, dat->dim, dat->type, acc, mapping);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, oppic_access acc, opp_mapping mapping) 
{
    if (dat == nullptr) 
    { 
        std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; 
        oppic_arg arg; 
        return arg; 
    }
    
    return oppic_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, acc, mapping);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, 
                                oppic_access acc, opp_mapping mapping) 
{
    oppic_arg arg;
    arg.index       = -1;
    arg.argtype     = OP_ARG_DAT;

    arg.dat         = dat;
    arg.map         = map;
    arg.dim         = dim;
    arg.idx         = idx;
    
    arg.size        = ((map == NULL) ? -1 : map->from->size + map->from->exec_size + map->from->nonexec_size);
    arg.data        = ((dat == NULL) ? NULL : dat->data);
    arg.data_d      = ((dat == NULL) ? NULL : dat->data_d);
    arg.map_data    = ((idx == -1 || idx == -2) ? NULL : map->map);
    arg.map_data_d  = ((idx == -1 || idx == -2) ? NULL : map->map_d);

    arg.type        = typ;     // Not used
    arg.acc         = acc;
    arg.opt         = 1;
    arg.sent        = 0;
    arg.mesh_mapping= mapping;
    
    return arg;
}

oppic_arg oppic_arg_dat_core(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, -1, NULL, acc, mapping);
}

// arg.map has the map, can change to mapping data map if required
oppic_arg oppic_arg_dat_core(oppic_map data_map, int idx, oppic_map map, oppic_access acc, 
                                opp_mapping mapping)
{
    oppic_arg arg;
    arg.argtype     = OP_ARG_MAP;

    arg.dat         = NULL;
    arg.map         = map;
    arg.dim         = data_map->dim;
    arg.idx         = idx;
    
    arg.size        = data_map->from->size + data_map->from->exec_size + data_map->from->nonexec_size;
    arg.data        = (char*)data_map->map;
    arg.data_d      = (char*)data_map->map_d;
    arg.map_data    = ((map == NULL) ? NULL : map->map);
    arg.map_data_d  = ((map == NULL) ? NULL : map->map_d);

    arg.type        = "int";
    arg.acc         = acc;
    arg.opt         = 1;

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
    arg.opt         = 1;

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
    arg.opt         = 1;

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
    arg.opt         = 1;
    
    return arg;
}

//****************************************
// oppic_decl_particle_set should make an oppic_set which changes size during the loop
oppic_set oppic_decl_particle_set_core(int size, char const *name, oppic_set cells_set) 
{   
    oppic_set set = oppic_decl_set_core(size, name);

    set->is_particle = true;
    set->cells_set = cells_set;

    return set;
}

// oppic_decl_particle_set should make an oppic_set which changes size during the loop
oppic_set oppic_decl_particle_set_core(char const *name, oppic_set cells_set) 
{    
    oppic_set set = oppic_decl_particle_set_core(0, name, cells_set);

    return set;
}

//****************************************
oppic_dat oppic_decl_particle_dat_core(oppic_set set, int dim, char const *type, int size, char *data, 
                                        char const *name, bool cell_index) 
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
bool oppic_increase_particle_count_core(oppic_set part_set, const int num_parts_to_insert)
{ 

    if (num_parts_to_insert <= 0) 
        return true;

    if (OP_DEBUG) 
        opp_printf("oppic_increase_particle_count_core", "set [%s] with size [%d]", 
            part_set->name, num_parts_to_insert);

    int new_part_set_size = part_set->size + num_parts_to_insert;

    // if the new particle set size is less or equal to set capacity, then just set new sizes instead of resizing
    if (part_set->set_capacity >= new_part_set_size)
    {
        if (OP_DEBUG) 
            opp_printf("opp_increase_particle_count_core", "set [%s] No need to realloc, new size[%d] set_capacity[%d]", 
                part_set->name, new_part_set_size, part_set->set_capacity);        
        
        part_set->size = new_part_set_size;
        part_set->diff = num_parts_to_insert;   
        return true;
    }

    // if the set needs resizing, then use alloc multiple to increase set capacity to reduce regular resizing of the set
    size_t new_part_set_capacity = part_set->size + num_parts_to_insert * OP_part_alloc_mult;
    bool return_flag = true;

    if (OP_DEBUG) //  || part_set->size != 0
        opp_printf("opp_increase_particle_count_core", "new_set_capacity %zu set_size %d num_dats_in_set %d", 
            new_part_set_capacity, part_set->size, part_set->particle_dats->size());

    // iterate over all the particle dats of that set and resize the arrays as necessary
    for (auto& dat : *(part_set->particle_dats))
    {
        if (dat->data == NULL) 
        {
            dat->data = (char *)malloc((size_t)(new_part_set_capacity * dat->size));
            // opp_printf("oppic_increase_particle_count_core", "malloc name %s %p size %d", 
            //     dat->name, dat->data, (new_part_set_capacity * dat->size));

            if (dat->data == nullptr)
            {
                opp_printf("opp_increase_particle_count_core", "Error... alloc of dat name %s failed (size %zu)", 
                    dat->name, (size_t)(new_part_set_capacity * dat->size));
                return_flag = false;
            }
        }
        else
        {
            // char* old = dat->data;
            dat->data = (char *)realloc(dat->data, (size_t)(new_part_set_capacity * dat->size));
            // opp_printf("oppic_increase_particle_count_core", "realloc %p name %s %p size %d", 
            //     old, dat->name, dat->data, (new_part_set_capacity * dat->size));

            if (dat->data == nullptr)
            {
                opp_printf("opp_increase_particle_count_core", "Error... realloc of dat name %s failed (size %zu)", 
                    dat->name, (size_t)(new_part_set_capacity * dat->size));
                return_flag = false;
            }
        }

        if (dat->is_cell_index && (dat->data != nullptr))
        {
            int* mesh_rel_array = (int *)dat->data;
            for (size_t i = (size_t)part_set->size; i < new_part_set_capacity; i++)
                mesh_rel_array[i] = MAX_CELL_INDEX;
        }
        // Note : The remainder in array from set size to set capacity will be garbage, except in mesh_relation_dat
    }
    
    part_set->size         = new_part_set_size;
    part_set->set_capacity = new_part_set_capacity;
    part_set->diff         = num_parts_to_insert;

    return return_flag;
}

//****************************************
// Not quite necessary - 
// If oppic_reset_num_particles_to_insert isn't called, set->diff will not reset and oppic_par_looppic_inject__ will 
// loop only on the injected particles; However it will reset upon oppic_increase_particle_count call
void oppic_reset_num_particles_to_insert_core(oppic_set set)
{
    set->diff = 0;
}

//****************************************
void oppic_init_particle_move_core(oppic_set set)
{
    if (OP_DEBUG) opp_printf("oppic_init_particle_move_core", "set [%s]", set->name);

    OPP_part_cells_set_size = set->cells_set->size;
    
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
        if (OP_DEBUG) 
            opp_printf("opp_mark_particle_to_move_core", "set [%s] particle_index [%d]", set->name, particle_index);

        set->particle_statuses[particle_index] = OPP_NEED_REMOVE;
        (set->particle_remove_count)++;
    }
    else if (move_status != (int)OPP_MOVE_DONE) 
    {
        std::cerr << "opp_mark_particle_to_move_core Failed to find the cell - Particle Index " << 
            particle_index << std::endl;
    }
}

#ifdef USE_OLD_FINALIZE
//****************************************
void oppic_finalize_particle_move_core(oppic_set set)
{
    if (OP_DEBUG) 
        opp_printf("oppic_finalize_particle_move_core", "set [%s] size[%d] with particle_remove_count [%d]", 
            set->name, set->size, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        // getting a backup of cell index since it will also be rearranged using a random OMP thread
        int *mesh_relation_data = (int *)malloc(set->set_capacity * set->mesh_relation_dat->size); 
        memcpy((char*)mesh_relation_data, set->mesh_relation_dat->data, 
                    set->set_capacity * set->mesh_relation_dat->size);

        for (int i = 0; i < (int)set->particle_dats->size(); i++)
        {
            oppic_dat current_oppic_dat = set->particle_dats->at(i);
            int removed_count = 0;
            int skip_count = 0;

            for (size_t j = 0; j < (size_t)set->size; j++)
            {
                if (set->particle_remove_count == (removed_count + skip_count))
                {
                    if (OP_DEBUG && i == 0) 
                        opp_printf("oppic_finalize_particle_move_core", 
                        "Required number already removed %d [%s] j=%d", 
                        set->particle_remove_count, set->name, j);
                    break;
                }

                if (mesh_relation_data[j] != MAX_CELL_INDEX) continue;

                char* dat_removed_ptr = (char *)(current_oppic_dat->data + (size_t)(j * current_oppic_dat->size));

                // BUG_FIX: (set->size - removed_count - 1) This index could marked to be removed, and if marked, 
                // then there could be an array index out of bounds access error in the future
                while ((set->size - removed_count - skip_count - 1 >= 0) && 
                    (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX))
                {
                    skip_count++;
                }
                if (j >= (size_t)(set->size - removed_count - skip_count - 1)) 
                {
                    if (OP_DEBUG && i == 0) opp_printf("oppic_finalize_particle_move_core", 
                        "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
                        j, (set->size - removed_count - skip_count - 1), set->name);
                    break;
                }
                if (set->size - removed_count - skip_count - 1 < 0)
                {
                    //if (OP_DEBUG && i == 0) 
                        opp_printf("oppic_finalize_particle_move_core", 
                        "(set->size - removed_count - skip_count - 1 < 0) %d %d %d [%s] j=%d %d", 
                        set->size, removed_count, skip_count, set->name, j, set->particle_remove_count);
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
            //     (size_t)(set->size - removed_count) * (size_t)current_oppic_dat->size);
        }

        free(mesh_relation_data);
    }
    else
    {
        if (OP_DEBUG) 
            opp_printf("oppic_finalize_particle_move_core", "Not processing dats since OP_auto_sort = TRUE");
    }

    set->size -= set->particle_remove_count;
}

#else

//****************************************
void oppic_finalize_particle_move_core(oppic_set set)
{
    if (OP_DEBUG) 
        opp_printf("oppic_finalize_particle_move_core", "set [%s] size[%d] with particle_remove_count [%d]", 
        set->name, set->size, set->particle_remove_count);

    // return if there are no particles to be removed
    if (set->particle_remove_count <= 0) 
        return;

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
        for (oppic_dat& dat : *(set->particle_dats))
        {
            for (const auto& x : swap_indices)
            {
                char* dat_removed_ptr = (char *)(dat->data + (x.first * dat->size));

                size_t offset_byte = x.second * dat->size;
                char* dat_to_replace_ptr = (char *)(dat->data + offset_byte);
                
                memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 
            }
        }
    }
    else
    {
        if (OP_DEBUG) 
            opp_printf("oppic_finalize_particle_move_core", "Not processing dats since OP_auto_sort = TRUE");
    }

    set->size -= set->particle_remove_count;
}
#endif

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

    if (OP_DEBUG) 
        printf("\toppic_remove_marked_particles_from_set set [%s] with size [%d]\n", 
            set->name, num_particles_to_remove);

    if (num_particles_to_remove <= 0) return;

    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {
        oppic_dat dat = set->particle_dats->at(i);
        int removed_count = 0;

        for (const int& j : idx_to_remove)
        {
            if (j < 0 || j > set->size) 
            {
                std::cerr << "oppic_remove_marked_particles_from_set set [" << set->name << 
                    "] has idx_to_remove [" << j << "]" << std::endl;
                continue;
            }

            char* dat_removed_ptr = (char *)(dat->data + (j * dat->size));
            // TODO : This will have a bug if the dat_to_replace is also to be removed! 
            // CHECK : oppic_finalize_particle_move_core and correct this if being used
            char* dat_to_replace_ptr = (char *)(dat->data + ((set->size - removed_count - 1) * dat->size));
            
            // Get the last element and replace the hole // Not the Optimum!!!
            // TODO : Can we make NULL data and handle it in sort?
            memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 

            removed_count++;
        }

        dat->data = (char *)realloc(dat->data, (size_t)(set->size - removed_count) * (size_t)dat->size);
    }

    set->size -= num_particles_to_remove;

    idx_to_remove.clear();
}

//****************************************
void oppic_particle_sort_core(oppic_set set)
{ 
    
    if (OP_DEBUG) printf("\toppic_particle_sort set [%s]\n", set->name);
    
    opp_profiler->start("PartSort");

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

    opp_profiler->end("PartSort");
}

//****************************************
void oppic_print_dat_to_txtfile_core(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{ 
    opp_profiler->start("PrintFile");

    const std::string file_name = std::string("files/") + file_name_prefix + "_" + file_name_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) 
    {
        opp_printf("oppic_print_dat_to_txtfile_core", "can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d -- %d %d\n", dat->set->size, dat->dim, dat->set->exec_size, dat->set->nonexec_size) < 0)
    {
        opp_printf("oppic_print_dat_to_txtfile_core", "error writing to %s\n", file_name.c_str());
        exit(2);
    }

    for (int i = 0; i < dat->set->size + dat->set->exec_size + dat->set->nonexec_size; i++) 
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

                if (fprintf(fp, " %2.25lE", ((double *)dat->data)[i * dat->dim + j]) < 0) 
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
                if (fprintf(fp, " %f", ((float *)dat->data)[i * dat->dim + j]) < 0) 
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
                if (fprintf(fp, " %d", ((int *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if ((strcmp(dat->type, "long") == 0) ||
                        (strcmp(dat->type, "long:soa") == 0)) 
            {
                if (fprintf(fp, " %ld", ((long *)dat->data)[i * dat->dim + j]) < 0) 
                {
                    printf("\toppic_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else 
            {
                printf("\toppic_print_dat_to_txtfile_core Unknown type %s, cannot be written to file %s\n", 
                    dat->type, file_name.c_str());
                exit(2);
            }
        }

        fprintf(fp, "\n");

        if (i+1 == dat->set->size) 
            fprintf(fp, "import_exec_below ****************************************\n");
        if (i+1 == dat->set->size + dat->set->exec_size) 
            fprintf(fp, "import_non_exec_below ****************************************\n");
    }
    
    fclose(fp);

    opp_profiler->end("PrintFile");
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
{ 
    opp_profiler->start("PrintFile");

    const std::string file_name = std::string("files/") + file_name_prefix + "_" + file_name_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) 
    {
        printf("\toppic_print_map_to_txtfile_core can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d -- %d %d\n", 
                map->from->size, map->dim, map->from->exec_size, map->from->nonexec_size) < 0) 
    {
        printf("\toppic_print_map_to_txtfile_core error writing to %s\n", file_name.c_str());
        exit(2);
    }

    for (int i = 0; i < map->from->size + map->from->exec_size + map->from->nonexec_size; i++) 
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

        if (i+1 == map->from->size) 
            fprintf(fp, "import_exec_below ****************************************\n");
        if (i+1 == map->from->size + map->from->exec_size) 
            fprintf(fp, "import_non_exec_below ****************************************\n");
    }
    
    fclose(fp);

    opp_profiler->end("PrintFile");
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
                if (fscanf(fp, " %lf %lf %lf %lf\n", 
                    &d_data[n * dim + 0], &d_data[n * dim + 1], &d_data[n * dim + 2], &d_data[n * dim + 3]) != 4) 
                        is_error = true;
                break;    
            case 16:
                if (fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                    &d_data[n * dim + 0], &d_data[n * dim + 1], &d_data[n * dim + 2], &d_data[n * dim + 3], 
                    &d_data[n * dim + 4], &d_data[n * dim + 5], &d_data[n * dim + 6], &d_data[n * dim + 7],
                    &d_data[n * dim + 8], &d_data[n * dim + 9], &d_data[n * dim + 10], &d_data[n * dim + 11], 
                    &d_data[n * dim + 12], &d_data[n * dim + 13], &d_data[n * dim + 14], &d_data[n * dim + 15]) != 16) 
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
                if (fscanf(fp, " %d %d %d\n", 
                    &i_data[n * dim + 0], &i_data[n * dim + 1], &i_data[n * dim + 2]) != 3) 
                        is_error = true;
                break;
            case 4:
                if (fscanf(fp, " %d %d %d %d\n", 
                    &i_data[n * dim + 0], &i_data[n * dim + 1], &i_data[n * dim + 2], &i_data[n * dim + 3]) != 4) 
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

//****************************************
bool opp_inc_part_count_with_distribution_core(oppic_set set, int num_parts_to_insert, oppic_dat part_dist)
{
    opp_profiler->start("IncPartCountWithDistribution");

    if (!oppic_increase_particle_count_core(set, num_parts_to_insert))
    {
        opp_profiler->end("IncPartCountWithDistribution");
        return false;
    }

    int* mesh_connectivity = (int *)set->mesh_relation_dat->data;
    int* distribution           = (int *)part_dist->data;

    int start = (set->size - set->diff);
    int j = 0;

    for (int i = 0; i < set->diff; i++)
    {
        if (i >= distribution[j]) j++; // check whether it is j or j-1    
        mesh_connectivity[start + i] = j;
    } 

    opp_profiler->end("IncPartCountWithDistribution");

    return true;
}

//****************************************
// Set Dirty::Device, making HOST data to be clean
void opp_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
            (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
            args[n].acc == OP_RW)) 
        {
            if (OP_DEBUG) opp_printf("opp_set_dirtybit", "Setting Dirty::Device| %s", args[n].dat->name);
            args[n].dat->dirty_hd = Dirty::Device;
            if (!args[n].dat->set->is_particle)
                args[n].dat->dirtybit = 1;
        }
    }
}

//****************************************
// Set Dirty::Host, making DEVICE data to be clean
void opp_set_dirtybit_cuda(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
            (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
            args[n].acc == OP_RW)) 
        {
            if (OP_DEBUG) opp_printf("opp_set_dirtybit_cuda", "Setting Dirty::Host| %s", args[n].dat->name);          
            args[n].dat->dirty_hd = Dirty::Host;
            if (!args[n].dat->set->is_particle)
                args[n].dat->dirtybit = 1;
        }
    }
}

//****************************************
void opp_set_dirtybit_grouped(int nargs, oppic_arg *args, DeviceType device)
{
    return (device == Device_CPU) ? opp_set_dirtybit(nargs, args) : opp_set_dirtybit_cuda(nargs, args);
}




//*******************************************************************************
//*******************************************************************************

#define BOUNDING_TOLERENCE 1e-12

//*******************************************************************************
BoundingBox::BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate, const std::shared_ptr<Comm> comm) {

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(dim, 0, comm);

    if (OP_DEBUG)
        opp_printf("Local BoundingBox [provided]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);

}

// For now, only 3D is implemented
//*******************************************************************************
BoundingBox::BoundingBox(const opp_dat node_pos_dat, int dim, const std::shared_ptr<Comm> comm) {
    
    if (dim != 3) {
        opp_abort(std::string("For now, only 3D BoundingBox is implemented"));
    }

    const double* node_pos_data = (const double*)node_pos_dat->data;
    const int node_count = node_pos_dat->set->size + node_pos_dat->set->exec_size + 
                                node_pos_dat->set->nonexec_size;;

    opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
    opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

    // make the bounding box even over the halo regions
    for (int i = 0; i < node_count; i++) {
        minCoordinate.x = std::min(node_pos_data[i * dim + 0], minCoordinate.x);
        minCoordinate.y = std::min(node_pos_data[i * dim + 1], minCoordinate.y);
        minCoordinate.z = std::min(node_pos_data[i * dim + 2], minCoordinate.z);
        maxCoordinate.x = std::max(node_pos_data[i * dim + 0], maxCoordinate.x);
        maxCoordinate.y = std::max(node_pos_data[i * dim + 1], maxCoordinate.y);
        maxCoordinate.z = std::max(node_pos_data[i * dim + 2], maxCoordinate.z);
    }

    if (node_count != 0) {
        minCoordinate.x -= BOUNDING_TOLERENCE;
        minCoordinate.y -= BOUNDING_TOLERENCE;
        minCoordinate.z -= BOUNDING_TOLERENCE;
        maxCoordinate.x += BOUNDING_TOLERENCE;
        maxCoordinate.y += BOUNDING_TOLERENCE;
        maxCoordinate.z += BOUNDING_TOLERENCE;
    }

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(dim, node_count, comm);

    if (OP_DEBUG)
        opp_printf("Local BoundingBox [computed]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);
}

//*******************************************************************************
void BoundingBox::generateGlobalBoundingBox(int dim, int count, const std::shared_ptr<Comm> comm) {

#ifdef USE_MPI

    const double* localMin = reinterpret_cast<const double*>(&(this->boundingBox[0]));
    double* globalMin = reinterpret_cast<double*>(&(this->globalBoundingBox[0]));    
    MPI_Allreduce(localMin, globalMin, dim, MPI_DOUBLE, MPI_MIN, comm->comm_parent);

    const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
    double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
    MPI_Allreduce(localMax, globalMax, dim, MPI_DOUBLE, MPI_MAX, comm->comm_parent);

    // This is to avoid min and max corrdinates to not have MAX_REAL and MIN_REAL when current rank has no work
    if (count == 0) {
        this->boundingBox[0].x = this->globalBoundingBox[0].x; 
        this->boundingBox[0].y = this->globalBoundingBox[0].y; 
        this->boundingBox[0].z = this->globalBoundingBox[0].z; 

        this->boundingBox[1].x = this->globalBoundingBox[0].x;
        this->boundingBox[1].y = this->globalBoundingBox[0].y;
        this->boundingBox[1].z = this->globalBoundingBox[0].z;
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("Global BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->globalBoundingBox[0].x, this->globalBoundingBox[0].y, this->globalBoundingBox[0].z, 
            this->globalBoundingBox[1].x, this->globalBoundingBox[1].y, this->globalBoundingBox[1].z);
#else
    this->globalBoundingBox = this->boundingBox;
#endif
}

//*******************************************************************************
BoundingBox::~BoundingBox() { };

//*******************************************************************************
const opp_point& BoundingBox::getLocalMin() const {

    return this->boundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getLocalMax() const {

    return this->boundingBox[1];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMin() const {

    return this->globalBoundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMax() const {
    
    return this->globalBoundingBox[1];
}

//*******************************************************************************
bool BoundingBox::isCoordinateInBoundingBox(const opp_point& point) { 

    if (this->boundingBox[0].x > point.x || this->boundingBox[1].x < point.x) 
        return false;
    else if (this->boundingBox[0].y > point.y || this->boundingBox[1].y < point.y) 
        return false;
    else if (this->boundingBox[0].z > point.z || this->boundingBox[1].z < point.z) 
        return false;
    
    return true;
}

//*******************************************************************************
bool BoundingBox::isCoordinateInGlobalBoundingBox(const opp_point& point) { 

    if (this->globalBoundingBox[0].x > point.x || this->globalBoundingBox[1].x < point.x) 
        return false;
    else if (this->globalBoundingBox[0].y > point.y || this->globalBoundingBox[1].y < point.y) 
        return false;
    else if (this->globalBoundingBox[0].z > point.z || this->globalBoundingBox[1].z < point.z) 
        return false;
    
    return true;
}


//*******************************************************************************
// This will contain mappings for halo indices too
GlobalToLocalCellIndexMapper::GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat) {

    // if (OP_DEBUG) 
    if (OPP_rank == 0)            
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping start");
        
    globalToLocalCellIndexMap.clear();

    const opp_set cells_set = global_cell_id_dat->set;
    int size_inc_halo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;

    for (int i = 0; i < size_inc_halo; i++) {

        int glbIndex = ((int*)global_cell_id_dat->data)[i];
        
        globalToLocalCellIndexMap.insert(std::make_pair(glbIndex, i));
    }
    
    // if (OP_DEBUG) 
    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping end");
}

//*******************************************************************************
GlobalToLocalCellIndexMapper::~GlobalToLocalCellIndexMapper() {
    
    globalToLocalCellIndexMap.clear();

    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "destroyGlobalToLocalCellIndexMapping");
}

//*******************************************************************************
int GlobalToLocalCellIndexMapper::map(const int globalIndex) {

    if (globalIndex == MAX_CELL_INDEX)
        return MAX_CELL_INDEX;

    auto it = globalToLocalCellIndexMap.find(globalIndex);
    if (it != globalToLocalCellIndexMap.end())
        return it->second;
    
    opp_printf("GlobalToLocalCellIndexMapper", "Error... local cell index not found for global index %d", globalIndex);
    return MAX_INT;
}


//*******************************************************************************
CellMapper::CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, const std::shared_ptr<Comm> comm) 
    : boundingBox(boundingBox), gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing), 
        minGlbCoordinate(boundingBox->getGlobalMin()), comm(comm) {
    
    const opp_point& minGblCoordinate = boundingBox->getGlobalMin();
    const opp_point& maxGblCoordinate = boundingBox->getGlobalMax();

    int ax = 0, ay = 0, az = 0;

    { // removed this and added below due to decimal point issues
        // this->globalGridDims.x = std::ceil((maxGblCoordinate.x - minGblCoordinate.x) * oneOverGridSpacing);
        // this->globalGridDims.y = std::ceil((maxGblCoordinate.y - minGblCoordinate.y) * oneOverGridSpacing);
        // this->globalGridDims.z = std::ceil((maxGblCoordinate.z - minGblCoordinate.z) * oneOverGridSpacing);   
    }
    {
        for (double z = minGblCoordinate.z; z < maxGblCoordinate.z; z += this->gridSpacing) az++;
        for (double y = minGblCoordinate.y; y < maxGblCoordinate.y; y += this->gridSpacing) ay++;
        for (double x = minGblCoordinate.x; x < maxGblCoordinate.x; x += this->gridSpacing) ax++; 
        this->globalGridDims.x = ax + 1;
        this->globalGridDims.y = ay + 1;
        this->globalGridDims.z = az + 1; 
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("CellMapper", "Global Grid Size - [%d %d %d] gridSpacing [%2.10lE]", 
            this->globalGridDims.x, this->globalGridDims.y, this->globalGridDims.z, this->gridSpacing); 
    
    const opp_point& minLocalCoordinate = boundingBox->getLocalMin();
    const opp_point& maxLocalCoordinate = boundingBox->getLocalMax();

    // Find the local ranks grid start indices
    ax = 0; ay = 0; az = 0;
    for (double z = minGblCoordinate.z; (z < minLocalCoordinate.z); z += this->gridSpacing) az++;
    for (double y = minGblCoordinate.y; (y < minLocalCoordinate.y); y += this->gridSpacing) ay++; 
    for (double x = minGblCoordinate.x; (x < minLocalCoordinate.x); x += this->gridSpacing) ax++; 
    this->localGridStart.x = ax == 0 ? ax : (ax - 1);
    this->localGridStart.y = ay == 0 ? ay : (ay - 1);
    this->localGridStart.z = az == 0 ? az : (az - 1);         

    // Find the local ranks grid end indices
    ax = 0; ay = 0; az = 0;
    for (double z = minGblCoordinate.z; (z <= maxLocalCoordinate.z); z += this->gridSpacing) az++; 
    for (double y = minGblCoordinate.y; (y <= maxLocalCoordinate.y); y += this->gridSpacing) ay++; 
    for (double x = minGblCoordinate.x; (x <= maxLocalCoordinate.x); x += this->gridSpacing) ax++; 
    this->localGridEnd.x = this->globalGridDims.x == ax ? ax : (ax + 1);
    this->localGridEnd.y = this->globalGridDims.y == ay ? ay : (ay + 1);
    this->localGridEnd.z = this->globalGridDims.z == az ? az : (az + 1);     

    if (OP_DEBUG)
        opp_printf("CellMapper", "Local Grid - Min[%d %d %d] Max[%d %d %d]", 
            this->localGridStart.x, this->localGridStart.y, this->localGridStart.z, 
            this->localGridEnd.x, this->localGridEnd.y, this->localGridEnd.z); 
}

//*******************************************************************************
CellMapper::~CellMapper() { 

#ifdef USE_MPI
    CHECK(MPI_Win_free(&this->win_structMeshToCellMapping))
    this->structMeshToCellMapping = nullptr;

    CHECK(MPI_Win_free(&this->win_structMeshToRankMapping))
    this->structMeshToRankMapping = nullptr;           
#else
    delete[] this->structMeshToCellMapping;
#endif
};

//*******************************************************************************
opp_point CellMapper::getCentroidOfBox(const opp_point& coordinate) { 

    opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);
    const opp_point& maxCoordinate = boundingBox->getGlobalMax();

    constexpr int DIM = 3; // Hardcoded for now
    switch (DIM) {
        case 1:
            ASSIGN_CENTROID_TO_DIM(x); 
            break;
        case 2:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            break;
        case 3:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            ASSIGN_CENTROID_TO_DIM(z);
            break;
        default:
            std::cerr << "Error getCentroidOfBox: Dimension invalid " << DIM << std::endl;
    }

    return centroid;
}

//*******************************************************************************
// Returns the global cell index
size_t CellMapper::findStructuredCellIndex(const opp_point& position) { 

    // Perform the calculations in higher precision (double)
    double xDiff = position.x - this->minGlbCoordinate.x;
    double yDiff = position.y - this->minGlbCoordinate.y;
    double zDiff = position.z - this->minGlbCoordinate.z;

    xDiff = xDiff * this->oneOverGridSpacing;
    yDiff = yDiff * this->oneOverGridSpacing;
    zDiff = zDiff * this->oneOverGridSpacing;

    // Round to the nearest integer to minimize rounding errors
    const int xIndex = static_cast<int>(xDiff);
    const int yIndex = static_cast<int>(yDiff);
    const int zIndex = static_cast<int>(zDiff);

    // Calculate the cell index mapping index
    size_t index = ((size_t)(xIndex) + (size_t)(yIndex * globalGridDims.x) + 
                (size_t)(zIndex * globalGridDims.x * globalGridDims.y));

    if (index >= globalGridSize) {
        // opp_printf("CellMapper", "Error index %zu generated is larger than globalGridSize %zu", 
        //     index, globalGridSize);
        return MAX_CELL_INDEX;
    }

    return index;
}

//*******************************************************************************
// Returns the global cell index
int CellMapper::findClosestCellIndex(const size_t& structCellIdx) { 
    
    if (OP_DEBUG) 
    {
        if (structCellIdx >= globalGridSize) {
            opp_printf("findClosestCellIndex", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                structCellIdx, globalGridSize);
            return MAX_CELL_INDEX;
        }
    }
    
    return this->structMeshToCellMapping[structCellIdx];
}

//*******************************************************************************
// Returns the rank of the cell
int CellMapper::findClosestCellRank(const size_t& structCellIdx) { 

#ifdef USE_MPI 
    if (OP_DEBUG) 
    {
        if (structCellIdx >= globalGridSize) {
            opp_printf("findClosestCellRank", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                structCellIdx, globalGridSize);
            return MAX_CELL_INDEX;
        }
    }

    return this->structMeshToRankMapping[structCellIdx];
#else
    return OPP_rank;
#endif
}

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) {

#ifdef USE_MPI
    waitBarrier();

    if (comm->rank_intra == 0) { // comm->size_intra > 1 && 

        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

        if (OP_DEBUG) 
            opp_printf("CellMapper", "reduceInterNodeMappings Start size_inter %d", comm->size_inter);

        std::string sendLog = "Send Counts: ", recvLog = "Recv Counts: ";

        int totalSendCount = 0, totalRecvCount = 0;
        std::vector<int> sendCounts(comm->size_inter);
        std::vector<int> sendDisplacements(comm->size_inter);
        std::vector<int> recvCounts(comm->size_inter);
        std::vector<int> recvDisplacements(comm->size_inter);

        for (int i = 0; i < comm->size_inter; ++i) {
            const int remainder = (int)(this->globalGridSize % comm->size_inter);
            sendCounts[i] = (this->globalGridSize / comm->size_inter) + (i < remainder ? 1 : 0);
            sendDisplacements[i] = totalSendCount;
            totalSendCount += sendCounts[i];

            if (OP_DEBUG) 
                sendLog += std::to_string(sendCounts[i]) + "|Disp:" + std::to_string(sendDisplacements[i]) + " ";
        }

        for (int i = 0; i < comm->size_inter; ++i) {
            recvCounts[i] = sendCounts[comm->rank_inter];
            recvDisplacements[i] = totalRecvCount;
            totalRecvCount += recvCounts[i];

            if (OP_DEBUG) 
                recvLog += std::to_string(recvCounts[i]) + "|Disp:" + std::to_string(recvDisplacements[i]) + " ";
        }

        if (OP_DEBUG) 
        {
            opp_printf("CellMapper", "reduceInterNodeMappings totalSendCount %d : %s", totalSendCount, sendLog.c_str());
            opp_printf("CellMapper", "reduceInterNodeMappings totalRecvCount %d : %s", totalRecvCount, recvLog.c_str());
        }

        std::vector<int> cellMappingsRecv(totalRecvCount, MAX_CELL_INDEX-1);
        std::vector<int> ranksRecv(totalRecvCount, MAX_CELL_INDEX-1);

        const int INT_COUNT_PER_MSG = 10000;

        std::vector<MPI_Request> sendRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalSendCount = sendCounts[rank];
            int sendCount = 0, alreadySentCount = 0;

            int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {
                
                if (i != blocks-1) 
                    sendCount = INT_COUNT_PER_MSG;
                else 
                    sendCount = totalSendCount - alreadySentCount;

                // opp_printf("SEND", "blocks %d|%d disp %d sendCounts %d to rank %d | %d %d %d %d %d %d %d %d %d %d", 
                //     i, blocks,
                //     sendDisplacements[rank] + i * INT_COUNT_PER_MSG, sendCount, inter_ranks[rank],
                //     sb[0], sb[1], sb[2], sb[3], sb[4], sb[5], sb[6], sb[7], sb[8], sb[9]);

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                        rank, (10000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                        rank, (20000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                alreadySentCount += sendCount;
            }
        }

        std::vector<MPI_Request> recvRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalRecvCount = recvCounts[rank];
            int recvCount = 0, alreadyRecvCount = 0;

            int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {

                if (i != blocks-1) 
                    recvCount = INT_COUNT_PER_MSG;
                else 
                    recvCount = totalRecvCount - alreadyRecvCount;

                // opp_printf("RECV", "blocks %d|%d disp %d recvCounts %d from rank %d", i, blocks,
                //     recvDisplacements[rank] + i * INT_COUNT_PER_MSG, recvCount, inter_ranks[rank]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(ranksRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                        rank, (10000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(cellMappingsRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                        rank, (20000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                alreadyRecvCount += recvCount;
            }
        }

        // opp_printf("DISPLACEMENTS", "send %d recv %d", sendDisplacements[comm->rank_inter], recvDisplacements[comm->rank_inter]);

        // Copy own data
        for (int i = 0; i < recvCounts[comm->rank_inter]; i++) {

            ranksRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i];
            cellMappingsRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i];
        }

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);

        // printStructuredMesh(std::string("RECV_MAPPING") + std::to_string(callID), cellMappingsRecv.data(), totalRecvCount);
        // printStructuredMesh(std::string("RECV_RANKS") + std::to_string(callID), ranksRecv.data(), totalRecvCount);
        // MPI_Barrier(comm->comm_inter);

        const size_t recvCount = recvCounts[0];
        for (size_t i = 0; i < recvCount; i++) { // reduce to get common mapping on all inter node ranks

            int cellIndex = MAX_CELL_INDEX;
            for (int r = 0; r < comm->size_inter; r++) {

                if (cellIndex > cellMappingsRecv[i + r * recvCount]) {

                    cellIndex = cellMappingsRecv[i + r * recvCount];
                    cellMappingsRecv[i] = cellIndex;
                    ranksRecv[i] = ranksRecv[i + r * recvCount];
                }
            }
        }

        MPI_Barrier(comm->comm_inter);

        // printStructuredMesh(std::string("MAPPING_COMP") + std::to_string(callID), ranksRecv.data(), recvCount);
        // MPI_Barrier(comm->comm_inter);

        sendRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalSendCount = recvCount;
            int sendCount = 0, alreadySentCount = 0;

            int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {
                
                if (i != blocks-1) 
                    sendCount = INT_COUNT_PER_MSG;
                else 
                    sendCount = totalSendCount - alreadySentCount;

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(cellMappingsRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (30000 + i), 
                        comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(ranksRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (40000 + i), 
                        comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                alreadySentCount += sendCount;
            }
        }

        // Since we are about to write to data mapped by windows, release the lock here
        CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
        CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));

        recvRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalRecvCount = sendCounts[rank];
            int recvCount2 = 0, alreadyRecvCount = 0;

            int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {

                if (i != blocks-1) 
                    recvCount2 = INT_COUNT_PER_MSG;
                else 
                    recvCount2 = totalRecvCount - alreadyRecvCount;

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                        rank, (30000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                        rank, (40000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                alreadyRecvCount += recvCount2;
            }
        }

        // Copy own data
        for (int i = 0; i < sendCounts[comm->rank_inter]; i++) {

            this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i] = ranksRecv[i];
            this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i] = cellMappingsRecv[i];
        }

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);
        
        if (OP_DEBUG) 
            opp_printf("CellMapper", "reduceInterNodeMappings END");
    }

    waitBarrier();
#endif
}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OP_DEBUG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat);

#ifdef USE_MPI
    if (comm->rank_intra == 0) {
        
        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));

        for (size_t i = 0; i < globalGridSize; i++) {
            
            if (this->structMeshToCellMapping[i] != MAX_CELL_INDEX)
                this->structMeshToCellMapping[i] = (-1 * this->structMeshToCellMapping[i]);
        }

        CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
    }

    MPI_Win_fence(0, this->win_structMeshToCellMapping); 
#endif

    for (size_t i = 0; i < globalGridSize; i++) {
        
        if (this->structMeshToRankMapping[i] == OPP_rank) {
            
            const int globalCID = (-1 * this->structMeshToCellMapping[i]);
            
            if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {
                
                const int localCID = globalToLocalCellIndexMapper.map(globalCID);
                
                if (localCID != MAX_CELL_INDEX) {
                    this->structMeshToCellMapping[i] = localCID;
                }
                else {
                    opp_printf("CellMapper", "Error at convertToLocalMappings : structMeshToCellMapping at index %d is invalid [gcid:%d]", 
                        i, this->structMeshToCellMapping[i]);
                }
            }
        }
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);

    if (comm->rank_intra == 0) {

        MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MAX, comm->comm_inter);
    }

    waitBarrier();
#endif
    if (OP_DEBUG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::enrichStructuredMesh(const int index, const int cell_index, const int rank) {

#ifdef USE_MPI
    MPI_Put(&cell_index, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToCellMapping);
    MPI_Put(&rank, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToRankMapping);
#else
    this->structMeshToCellMapping[index] = cell_index;
#endif
}

//*******************************************************************************
void CellMapper::printStructuredMesh(const std::string msg, int *array, size_t size, bool printToFile) {
    
    // if (!OP_DEBUG)
    //     return;

    if (!printToFile) 
    {
        opp_printf("structMeshToCellMapping", "%s - size=%zu", msg.c_str(), size);

        for (size_t i = 0; i < size; i++) {
            printf("%zu|%d ", i, array[i]);
            if (i % 50 == 0) printf("\n");
        }   
        printf("\n");
    }
    else
    {
        const std::string file_name = std::string("files/struct_com") + std::to_string(OPP_comm_size) + "_r" + 
                                        std::to_string(OPP_rank) + "_" + msg; 

        std::ofstream outFile(file_name);

        if (!outFile.is_open()) {
            opp_printf("printStructuredMesh", "can't open file %s\n", file_name.c_str());
            opp_abort("printStructuredMesh - can't open file");
        }

        for (size_t i = 0; i < size; i++) {

            outFile << i << "|";

            int value = array[i];
            if (value != MAX_CELL_INDEX) {
                outFile << value << " ";
            }
            else {
                outFile << "X ";
            }

            if ((i+1) % 50 == 0) 
                outFile << "\n";
        }   
        outFile << "\n";
        outFile.close();
    }
}

//*******************************************************************************
void CellMapper::createStructMeshMappingArrays() {

    globalGridSize = (size_t)(globalGridDims.x * globalGridDims.y * globalGridDims.z);

#ifdef USE_MPI

    MPI_Barrier(MPI_COMM_WORLD);

    // One per shared memory (node)
    
    // create CELL INDEX array
    {
        const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

        CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                    (void *)&this->structMeshToCellMapping, &this->win_structMeshToCellMapping))

        MPI_Aint allocatedSize = 0;
        int disp = 0;
        CHECK(MPI_Win_shared_query(this->win_structMeshToCellMapping, 0, &allocatedSize, &disp,
                                    (void *)&this->structMeshToCellMapping))

        if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
            opp_abort(std::string("Pointer to incorrect size in MPI structMeshToCellMapping"));
        }
        if (disp != sizeof(int)) {
            opp_abort(std::string("Invalid displacement unit in MPI structMeshToCellMapping"));
        }

        if (comm->rank_intra == 0) {
            for (size_t i = 0; i < globalGridSize; i++)
                this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
        }

        MPI_Win_fence(0, this->win_structMeshToCellMapping);
    }

    // create RANK array
    {
        const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

        CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                    (void *)&this->structMeshToRankMapping, &this->win_structMeshToRankMapping))

        MPI_Aint allocatedSize = 0;
        int disp = 0;
        CHECK(MPI_Win_shared_query(this->win_structMeshToRankMapping, 0, &allocatedSize, &disp,
                                    (void *)&this->structMeshToRankMapping))

        if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
            opp_abort(std::string("Pointer to incorrect size in MPI structMeshToRankMapping"));
        }
        if (disp != sizeof(int)) {
            opp_abort(std::string("Invalid displacement unit in MPI structMeshToRankMapping"));
        }

        if (comm->rank_intra == 0) {
            for (size_t i = 0; i < globalGridSize; i++)
                this->structMeshToRankMapping[i] = MAX_CELL_INDEX;
        }

        MPI_Win_fence(0, this->win_structMeshToRankMapping);
    }
#else
    this->structMeshToCellMapping = new int[globalGridSize];
    for (size_t i = 0; i < globalGridSize; i++)
        this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
#endif
}

//*******************************************************************************
void CellMapper::waitBarrier() {

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, this->win_structMeshToCellMapping); 
    MPI_Win_fence(0, this->win_structMeshToRankMapping); 
#endif
}

//*******************************************************************************
void CellMapper::lockWindows() {

#ifdef USE_MPI
    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));
#endif
}


void CellMapper::unlockWindows() {

#ifdef USE_MPI
    CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
#endif
}

