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

#include "oppic.h"

std::list<op_set> op_sets;
std::list<op_map> op_maps;
std::list<op_dat> op_dats;

//****************************************
void op_init() {
    op_sets.clear();
    op_maps.clear();
    op_dats.clear();
}

//****************************************
void op_exit() {
    for (auto a : op_maps) {
        free(a->map);
        free(a);
    }
    op_maps.clear();

    for (auto a : op_dats) {
        free(a->data);
        for (int thr = 0; thr < (int)a->thread_data.size(); thr++) free(a->thread_data[thr]);
        free(a);
    }
    op_dats.clear();

    for (auto a : op_sets) {
        free(a);
    }
    op_sets.clear();

    #ifdef USE_TRACE
        std::string trace_file_name = std::string("trace/trace_") + getTimeStr() + ".csv";
        trace::current.write_profile(trace_file_name);
    #endif
}

//****************************************
static char *copy_str(char const *src) {
    const size_t len = strlen(src) + 1;
  char *dest = (char *)malloc(len * sizeof(char));
  
  return strncpy(dest, src, len);
}

//****************************************
op_set op_decl_set(int size, char const *name) {
    
    op_set set = (op_set)malloc(sizeof(op_set_core));
    set->size = size;
    set->name = copy_str(name);
    
    op_sets.push_back(set);
    return set;
}

//****************************************
op_map op_decl_map(op_set from, op_set to, int dim, int *imap, char const *name) {

    op_map map = (op_map)malloc(sizeof(op_map_core));
    map->from = from;
    map->to = to;
    map->dim = dim;
    // map->map = imap;
    map->map = (int *)malloc((size_t)from->size * (size_t)dim * sizeof(int)); // use op2's copy instead of imap;
    memcpy(map->map, imap, sizeof(int) * from->size * dim);
    map->name = copy_str(name);

    op_maps.push_back(map);
    return map;
}

//****************************************
op_dat op_decl_dat(op_set set, int dim, char const *type, int size, char *data, char const *name) {

    op_dat dat = (op_dat)malloc(sizeof(op_dat_core));
    dat->set = set;
    dat->dim = dim;
    dat->data = (char *)malloc((size_t)dim * (size_t)size * (size_t)(set->size) * sizeof(char));
    memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    dat->name = copy_str(name);
    dat->type = copy_str(type);
    dat->size = dim * size;

    op_dats.push_back(dat);
    return dat;
}

//****************************************
op_arg op_arg_dat(op_dat dat, int idx, op_map map, op_access acc, bool map_with_cell_index) {
    if (dat == nullptr) { std::cerr << "dat is NULL at op_arg_dat" << std::endl; op_arg arg; return arg; }
    return op_arg_dat(dat, idx, map, dat->dim, dat->type, acc, map_with_cell_index);
}

op_arg op_arg_dat(op_dat dat, op_access acc, bool map_with_cell_index) {
    if (dat == nullptr) { std::cerr << "dat is NULL at op_arg_dat" << std::endl; op_arg arg; return arg; }
    return op_arg_dat(dat, -1, NULL, dat->dim, dat->type, acc, map_with_cell_index);
}

op_arg op_arg_dat(op_dat dat, int idx, op_map map, int dim, const char *typ, op_access acc, bool map_with_cell_index) {

    op_arg arg;
    arg.argtype = OP_ARG_DAT;

    arg.dat = dat;
    arg.map = map;
    arg.dim = dim;
    arg.idx = idx;
    
    arg.size = ((dat == NULL) ? 0 : dat->size);
    arg.data = ((dat == NULL) ? 0 : dat->data);
    arg.map_data = ((idx == -1 || idx == -2) ? NULL : map->map);
    arg.type = typ;
    arg.acc = acc;

    return arg;
}

op_arg op_arg_dat(op_map map, op_access acc, bool map_with_cell_index)
{
    op_arg arg;
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
op_arg op_arg_gbl(double *data, int dim, char const *typ, op_access acc)
{
    op_arg arg;
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

op_arg op_arg_gbl(const bool *data, int dim, char const *typ, op_access acc)
{
    op_arg arg;
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


#ifdef OP_PARTICLES

//****************************************
op_set op_decl_particle_set(int size, char const *name, op_set cells_set) { // op_decl_particle_set should make and op_set which changes size during the loop
    
    op_set set = op_decl_set(size, name);
    set->is_particle = true;
    set->particle_dats.clear();
    set->cells_set = cells_set;
    return set;
}

//****************************************
op_set op_decl_particle_set(char const *name, op_set cells_set) { // op_decl_particle_set should make and op_set which changes size during the loop
    
    op_set set = op_decl_particle_set(0, name, cells_set);
    return set;
}

//****************************************
op_dat op_decl_particle_dat(op_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index) {

    op_dat dat = (op_dat)malloc(sizeof(op_dat_core));
    dat->set = set;
    dat->dim = dim;
    if (set->size > 0)
    {
        dat->data = (char *)malloc((size_t)dim * (size_t)size * (size_t)(set->size) * sizeof(char));
        memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    }

    dat->name = copy_str(name);
    dat->type = copy_str(type);
    dat->size = dim * size;
    dat->is_cell_index = cell_index;

    if (cell_index) set->cell_index_dat = dat;
    set->particle_dats.push_back(dat);
    op_dats.push_back(dat);

    return dat;
}

//****************************************
void op_increase_particle_count(op_set particles_set, const int num_particles_to_insert)
{
    if (num_particles_to_insert <= 0) return;

    if (OP_DEBUG) printf("op_increase_particle_count set [%s] with size [%d]\n", particles_set->name, num_particles_to_insert);

    particles_set->diff = num_particles_to_insert;
    particles_set->size += num_particles_to_insert;

    for (auto& current_op_dat : particles_set->particle_dats)
    {
        if (current_op_dat->data == nullptr) 
        {
            current_op_dat->data = (char *)malloc((size_t)(particles_set->size * current_op_dat->size));
        }
        else
        {
            current_op_dat->data = (char *)realloc(current_op_dat->data, (size_t)(particles_set->size * current_op_dat->size));
        }
    }  
}

//****************************************
// Not quite necessary - 
// If op_reset_num_particles_to_insert isn't called, set->diff will not reset and op_par_loop_inject__ will loop only on the injected particles
// However it will reset upon op_increase_particle_count call
void op_reset_num_particles_to_insert(op_set set)
{
    set->diff = 0;
}

//****************************************
void op_mark_particle_to_remove(op_set set, int particle_index)
{
    set->indexes_to_remove.push_back(particle_index);
}

//****************************************
void op_remove_marked_particles_from_set(op_set set)
{
    op_remove_marked_particles_from_set(set, set->indexes_to_remove);
}

void op_remove_marked_particles_from_set(op_set set, std::vector<int>& idx_to_remove)
{
    int num_particles_to_remove = idx_to_remove.size();

    std::sort(idx_to_remove.begin(), idx_to_remove.end()); // WE DONT NEED THIS IF NOT OMP OR CUDA

    if (OP_DEBUG) printf("op_remove_marked_particles_from_set set [%s] with size [%d]\n", set->name, num_particles_to_remove);

    if (num_particles_to_remove <= 0) return;

    for (int i = 0; i < (int)set->particle_dats.size(); i++)
    {
        op_dat current_op_dat = set->particle_dats[i];
        int removed_count = 0;

        for (const int& j : idx_to_remove)
        {
            if (j < 0 || j > set->size) 
            {
                std::cerr << "op_remove_marked_particles_from_set set [" << set->name << "] has idx_to_remove [" << j << "]" << std::endl;
                continue;
            }

            char* dat_removed_ptr = (char *)(current_op_dat->data + (j * current_op_dat->size));
            char* dat_to_replace_ptr = (char *)(current_op_dat->data + ((set->size - removed_count - 1) * current_op_dat->size));
            
            memcpy(dat_removed_ptr, dat_to_replace_ptr, current_op_dat->size); // Get the last element and replace the hole // Do we need to sort by cell index?? GOOD QUESTION!

            removed_count++;
        }

        current_op_dat->data = (char *)realloc(current_op_dat->data, (size_t)(set->size - removed_count) * (size_t)current_op_dat->size);
    }

    set->size -= num_particles_to_remove;

    idx_to_remove.clear();
}

//**************************************** // UNUSED
void op_insert_to_particle_dat(op_dat dat, const char* data) 
{
    char* data_copy = (char*)malloc((size_t)dat->dim * (size_t)dat->size * sizeof(char));
    memcpy(data_copy, data, (size_t)dat->dim * (size_t)dat->size * sizeof(char));
    dat->data_to_insert.push_back(data_copy);
}

//**************************************** // UNUSED
void op_insert_particles_to_set(op_set set) {

    int num_particles_to_insert = set->cell_index_dat->data_to_insert.size();

    if (num_particles_to_insert <= 0) return;

    bool is_valid = true;
    
    for (int i = 0; i < (int)set->particle_dats.size(); i++)
    {
        op_dat current_op_dat = set->particle_dats[i];

        if (num_particles_to_insert != (int)current_op_dat->data_to_insert.size())
        {
            is_valid = false;
            std::cerr << "op_insert_particles_to_set called with different inset dat vector sizes\n" << std::endl;
        }

        if (current_op_dat->data == nullptr) 
        {
            current_op_dat->data = (char *)malloc((size_t)num_particles_to_insert * (size_t)current_op_dat->size);
        }
        else
        {
            current_op_dat->data = (char *)realloc(current_op_dat->data, (size_t)(set->size + num_particles_to_insert) * (size_t)current_op_dat->size);
        }
    }

    if (OP_DEBUG) printf("op_insert_particles_to_set Allocating set [%s] with size [%d]\n", "PARTICLES", num_particles_to_insert);

    if (!is_valid) return;
                     
    for (int i = 0; i < (int)set->particle_dats.size(); i++)
    {
        op_dat current_op_dat = set->particle_dats[i];

        for (int j = 0; j < num_particles_to_insert; j++)
        {
            char* data_ptr = current_op_dat->data_to_insert[j];
            char* op_dat_ptr = (char *)(current_op_dat->data + ((set->size + j) * current_op_dat->size));
            
            memcpy(op_dat_ptr, data_ptr, current_op_dat->size);

            free(data_ptr);
        }

        current_op_dat->data_to_insert.clear();
    }

    set->size += num_particles_to_insert;
}

//****************************************
template <typename T>
std::vector<size_t> sort_indexes(const int* v, int size) 
{ TRACE_ME;

    std::vector<size_t> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
        [v](size_t i1, size_t i2) { return v[i1] < v[i2];});

    return idx;
}

//****************************************
void op_particle_sort(op_set set)
{ TRACE_ME;
    
    if (OP_DEBUG) printf("op_particle_sort set [%s]\n", set->name);
    
    int* cell_index_data = (int*)set->cell_index_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes<int>(cell_index_data, set->size);

    #pragma omp parallel for
    for (int i = 0; i < set->particle_dats.size(); i++)
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

            // for (auto& a : map) printf("%d -> start %d end %d |\t", a.first, a.second.start, a.second.end);
        }
    }
}

#endif

//*************************************************************************************************

#ifdef USE_TRACE
    int trace::enabled = 1;
    Trace trace::current = Trace("__TRACE_BASE__");
#endif

std::string getTimeStr()
{
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y_%m_%d__%H_%M_%S", std::localtime(&now));
    return s;
}