#include <oppic_lib_core.h>


#ifdef USE_TRACE
    int trace::enabled = 1;
    Trace trace::current = Trace("__TRACE_BASE__");
#endif


//****************************************
std::vector<oppic_set> oppic_sets;
std::vector<oppic_map> oppic_maps;
std::vector<oppic_dat> oppic_dats;

//****************************************
void oppic_init_core(int argc, char **argv, int diags) 
{
    oppic_sets.clear();
    oppic_maps.clear();
    oppic_dats.clear();
}

//****************************************
void oppic_exit_core() 
{
    for (auto a : oppic_maps) {
        free(a->map);
        free(a);
    }
    oppic_maps.clear();

    for (auto a : oppic_dats) {
        free(a->data);
        for (int thr = 0; thr < (int)a->thread_data.size(); thr++) { free(a->thread_data[thr]); }
        free(a);
    }
    oppic_dats.clear();

    for (auto a : oppic_sets) {
        free(a);
    }
    oppic_sets.clear();

    #ifdef USE_TRACE
        std::string trace_file_name = std::string("trace/trace_") + getTimeStr() + ".csv";
        trace::current.write_profile(trace_file_name);
    #endif
}

//****************************************
oppic_set oppic_decl_set_core(int size, char const *name) 
{
    oppic_set set = (oppic_set)malloc(sizeof(oppic_set_core));
    set->size = size;
    set->name = copy_str(name);
    
    oppic_sets.push_back(set);
    return set;
}

//****************************************
oppic_map oppic_decl_map_core(oppic_set from, oppic_set to, int dim, int *imap, char const *name) 
{
    oppic_map map = (oppic_map)malloc(sizeof(oppic_map_core));
    map->from = from;
    map->to = to;
    map->dim = dim;
    // map->map = imap;
    map->map = (int *)malloc((size_t)from->size * (size_t)dim * sizeof(int)); // use op2's copy instead of imap;
    memcpy(map->map, imap, sizeof(int) * from->size * dim);
    map->name = copy_str(name);

    oppic_maps.push_back(map);
    return map;
}

//****************************************
oppic_dat oppic_decl_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name) 
{
    oppic_dat dat = (oppic_dat)malloc(sizeof(oppic_dat_core));
    dat->set = set;
    dat->dim = dim;
    dat->data = (char *)malloc((size_t)dim * (size_t)size * (size_t)(set->size) * sizeof(char));
    memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    dat->name = copy_str(name);
    dat->type = copy_str(type);
    dat->size = dim * size;
    dat->thread_data.clear();

    oppic_dats.push_back(dat);
    return dat;
}

//****************************************
oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index) 
{
    if (dat == nullptr) { std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; oppic_arg arg; return arg; }
    
    return oppic_arg_dat_core(dat, idx, map, dat->dim, dat->type, acc, map_with_cell_index);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, oppic_access acc, bool map_with_cell_index) 
{
    if (dat == nullptr) { std::cerr << "dat is NULL at oppic_arg_dat" << std::endl; oppic_arg arg; return arg; }
    
    return oppic_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, acc, map_with_cell_index);
}

oppic_arg oppic_arg_dat_core(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index) 
{
    oppic_arg arg;
    arg.argtype = OP_ARG_DAT;

    arg.dat = dat;
    arg.map = map;
    arg.dim = dim;
    arg.idx = idx;
    
    arg.size = ((dat == NULL) ? 0 : dat->size);
    arg.data = ((dat == NULL) ? nullptr : dat->data);
    arg.map_data = ((idx == -1 || idx == -2) ? NULL : map->map);
    arg.type = typ;
    arg.acc = acc;

    return arg;
}

oppic_arg oppic_arg_dat_core(oppic_map map, oppic_access acc, bool map_with_cell_index)
{
    oppic_arg arg;
    arg.argtype = OP_ARG_DAT;

    arg.dat = NULL;
    arg.map = map;
    arg.dim = map->dim;
    arg.idx = -1;
    
    arg.size = map->from->size;
    arg.data = (char*)map->map;

    arg.type = "int";
    arg.acc = acc;

    return arg;
}

//****************************************
oppic_arg oppic_arg_gbl_core(double *data, int dim, char const *typ, oppic_access acc)
{
    oppic_arg arg;
    arg.argtype = OP_ARG_GBL;

    arg.dat = NULL;
    arg.map = NULL;
    arg.dim = dim;
    arg.idx = -1;
    arg.size = dim * sizeof(double);
    arg.data = (char*)data;
    arg.map_data = NULL;
    arg.type = typ;
    arg.acc = acc;
    
    return arg;
}

oppic_arg oppic_arg_gbl_core(const bool *data, int dim, char const *typ, oppic_access acc)
{
    oppic_arg arg;
    arg.argtype = OP_ARG_GBL;

    arg.dat = NULL;
    arg.map = NULL;
    arg.dim = dim;
    arg.idx = -1;
    arg.size = dim * sizeof(bool);
    arg.data = (char*)data;
    arg.map_data = NULL;
    arg.type = typ;
    arg.acc = acc;
    
    return arg;
}

//****************************************
oppic_set oppic_decl_particle_set_core(int size, char const *name, oppic_set cells_set) // oppic_decl_particle_set should make and oppic_set which changes size during the loop
{     
    oppic_set set = oppic_decl_set_core(size, name);

    set->is_particle = true;
    set->particle_dats.clear();
    set->cells_set = cells_set;
    set->diff = 0;
    set->cell_index_dat = nullptr;
    set->cell_index_v_part_index_map.clear();

    return set;
}

//****************************************
oppic_set oppic_decl_particle_set_core(char const *name, oppic_set cells_set) // oppic_decl_particle_set should make and oppic_set which changes size during the loop
{    
    oppic_set set = oppic_decl_particle_set_core(0, name, cells_set);

    return set;
}

//****************************************
oppic_dat oppic_decl_particle_dat_core(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index) {

    oppic_dat dat = (oppic_dat)malloc(sizeof(oppic_dat_core));
    dat->set = set;
    dat->dim = dim;
    if (set->size > 0)
    {
        dat->data = (char *)malloc((size_t)dim * (size_t)size * (size_t)(set->size) * sizeof(char));
        memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    }
    else
    {
        dat->data = nullptr;
    }

    dat->name = copy_str(name);
    dat->type = copy_str(type);
    dat->size = dim * size;
    dat->is_cell_index = cell_index;
    dat->thread_data.clear();

    if (cell_index) set->cell_index_dat = dat;
    set->particle_dats.push_back(dat);
    oppic_dats.push_back(dat);

    return dat;
}

//****************************************
void oppic_increase_particle_count_core(oppic_set particles_set, const int num_particles_to_insert)
{
    if (num_particles_to_insert <= 0) return;

    if (OP_DEBUG) printf("oppic_increase_particle_count set [%s] with size [%d]\n", particles_set->name, num_particles_to_insert);

    particles_set->diff = num_particles_to_insert;
    particles_set->size += num_particles_to_insert;

    for (auto& current_oppic_dat : particles_set->particle_dats)
    {
        if (current_oppic_dat->data == nullptr) 
        {
            current_oppic_dat->data = (char *)malloc((size_t)(particles_set->size * current_oppic_dat->size));
        }
        else
        {
            current_oppic_dat->data = (char *)realloc(current_oppic_dat->data, (size_t)(particles_set->size * current_oppic_dat->size));
        }
    }  
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
void oppic_mark_particle_to_remove_core(oppic_set set, int particle_index)
{
    set->indexes_to_remove.push_back(particle_index);
}

//****************************************
void oppic_remove_marked_particles_from_set_core(oppic_set set)
{
    oppic_remove_marked_particles_from_set_core(set, set->indexes_to_remove);
}

void oppic_remove_marked_particles_from_set_core(oppic_set set, std::vector<int>& idx_to_remove)
{
    int num_particles_to_remove = idx_to_remove.size();

    std::sort(idx_to_remove.begin(), idx_to_remove.end());

    if (OP_DEBUG) printf("oppic_remove_marked_particles_from_set set [%s] with size [%d]\n", set->name, num_particles_to_remove);

    if (num_particles_to_remove <= 0) return;

    for (int i = 0; i < (int)set->particle_dats.size(); i++)
    {
        oppic_dat current_oppic_dat = set->particle_dats[i];
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
{ 
    
    if (OP_DEBUG) printf("oppic_particle_sort set [%s]\n", set->name);
    
    int* cell_index_data = (int*)set->cell_index_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(cell_index_data, set->size);

    for (int i = 0; i < (int)set->particle_dats.size(); i++)
    {    
        auto& dat = set->particle_dats[i];
        char *new_data = (char *)malloc(set->size * dat->size);
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
            std::map<int, part_index>& map = set->cell_index_v_part_index_map;
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