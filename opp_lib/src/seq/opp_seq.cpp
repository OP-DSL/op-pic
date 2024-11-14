
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

#include <opp_seq.h>

opp_move_var move_var;

char opp_move_status_flag = OPPX_MOVE_DONE;
bool opp_move_hop_iter_one_flag = true;
OPP_INT* opp_p2c = nullptr;
OPP_INT* opp_c2c = nullptr;

//****************************************
void opp_init(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULLPTR, "opp::PetscSEQ");
#endif

    std::string log = "Running on SEQ";
    opp_printf("OP-PIC", "%s", log.c_str());
    opp_printf("OP-PIC", "---------------------------------------------");

    opp_init_core(argc, argv);
    opp_params->write(std::cout);
}

//****************************************
void opp_exit()
{

#ifdef USE_PETSC
    PetscFinalize();
#endif

    opp_exit_core();
}

//****************************************
void opp_abort(std::string s)
{
    opp_printf("opp_abort", "%s", s.c_str());
    exit(-1);
}

//****************************************
opp_set opp_decl_set(int size, char const *name)
{
    return opp_decl_set_core(size, name);
}

//****************************************
opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name)
{
    return opp_decl_map_core(from, to, dim, imap, name);
}

//****************************************
opp_dat opp_decl_mesh_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return opp_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);
}

//****************************************
opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    if (set->is_particle) 
        return opp_decl_part_dat(set, dim, dtype, data, name, cell_index);
    else
        return opp_decl_mesh_dat(set, dim, dtype, data, name);
}

//****************************************
opp_map opp_decl_map_txt(opp_set from, opp_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)opp_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    opp_map map = opp_decl_map(from, to, dim, map_data, name);

    opp_host_free(map_data);

    return map;
}

//****************************************
opp_dat opp_decl_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)opp_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    opp_dat dat = opp_decl_dat_core(set, dim, type.c_str(), size, dat_data, name);

    opp_host_free(dat_data);

    return dat;
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat1");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat2");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat3");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat4");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat5");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat6");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, nullptr, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_map data_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat7");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
                data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat8");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat9");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, idx, map, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat10");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, idx, map, p2c_map->p2c_dat, acc);
}

//****************************************
// template <class T> opp_arg opp_arg_gbl(T *data, int dim, char const *typ, opp_access acc);
opp_arg opp_arg_gbl(double *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}
opp_arg opp_arg_gbl(int *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}
opp_arg opp_arg_gbl(const bool *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set)
{
    return opp_decl_particle_set_core(size, name, cells_set);
}
opp_set opp_decl_particle_set(char const *name, opp_set cells_set)
{
    return opp_decl_particle_set_core(name, cells_set);
}

//****************************************
opp_dat opp_decl_part_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return opp_decl_particle_dat_core(set, dim, type.c_str(), size, (char*)data, name, cell_index);
}

//****************************************
opp_dat opp_decl_particle_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)opp_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    opp_dat dat = opp_decl_particle_dat_core(set, dim, type.c_str(), size, dat_data, name, cell_index);

    opp_host_free(dat_data);

    return dat;
}

//****************************************
void opp_increase_particle_count(opp_set particles_set, const int num_particles_to_insert)
{
    if (!opp_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_printf("opp_increase_particle_count", "Error: opp_increase_particle_count_core failed for particle set [%s]", particles_set->name);
        exit(-1);        
    }
}

//****************************************
void opp_inc_part_count_with_distribution(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist, bool calc_new)
{
    if (OPP_DBG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        opp_abort("opp_inc_part_count_with_distribution_core");        
    }
}

//****************************************
void opp_reset_num_particles_to_insert(opp_set set)
{
    opp_reset_num_particles_to_insert_core(set);
}

//****************************************
void opp_mark_particle_to_remove(opp_set set, int particle_index)
{
    opp_mark_particle_to_remove_core(set, particle_index);
}

//****************************************
void opp_remove_marked_particles_from_set(opp_set set)
{
    opp_remove_marked_particles_from_set_core(set);
}
void opp_remove_marked_particles_from_set(opp_set set, std::vector<int>& idx_to_remove)
{
    opp_remove_marked_particles_from_set_core(set, idx_to_remove);
}

//****************************************
void opp_particle_sort(opp_set set)
{ 
    opp_particle_sort_core(set);
}

//****************************************
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_s";
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    opp_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

//****************************************
void opp_dump_dat(opp_dat dat)
{
    opp_dump_dat_core(dat);
}

//****************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{ 
    opp_init_particle_move_core(set);

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
}

//****************************************
void opp_mark_particle_to_move(opp_set set, int particle_index, int move_status)
{
    // Can fill the holes here, since there could be removed particles
    opp_mark_particle_to_move_core(set, particle_index, move_status);
}

//****************************************
bool opp_finalize_particle_move(opp_set set)
{ 
    opp_finalize_particle_move_core(set);

    if (OPP_auto_sort == 1)
    {
        if (OPP_DBG) printf("\topp_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        opp_particle_sort(set);
    }

    // in seq, no need to iterate the communication loop
    return false;
}

//****************************************
void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset)
{
    int set_size = dat->set->size;

    for (int i = 0; i < set_size; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
int opp_mpi_halo_exchanges_grouped(opp_set set, int nargs, opp_arg *args, DeviceType device) 
{
    return set->size;
}

//****************************************
void opp_mpi_force_halo_update_if_dirty(opp_set set, std::vector<opp_dat> dats, DeviceType device) 
{
    // Nothing to execute here
}

//****************************************
void opp_mpi_halo_wait_all(int nargs, opp_arg *args)
{
    // Nothing to execute here
}

//****************************************
void opp_init_double_indirect_reductions(int nargs, opp_arg *args)
{
    // Nothing to execute here
}

//****************************************
void opp_exchange_double_indirect_reductions(int nargs, opp_arg *args)
{
    // Nothing to execute here
}


void opp_complete_double_indirect_reductions(int nargs, opp_arg *args)
{
    // Nothing to execute here
}

//****************************************
bool opp_part_check_status(opp_move_var& m, int map0idx, opp_set set, 
    int particle_index, int& remove_count, int thread)
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        remove_count += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }

    return true;
}

//****************************************
opp_move_var opp_get_move_var(int thread)
{
    move_var.move_status = OPP_MOVE_DONE;
    move_var.iteration_one = true;

    // no perf improvement by using a buffered move var, could create a new here instead
    
    return move_var;
}

//*******************************************************************************
// Copy a dat from host to device
void opp_upload_dat(opp_dat dat) {}

//*******************************************************************************
// Copy a dat from device to host
void opp_download_dat(opp_dat dat) {}

//*******************************************************************************
// Copy all dats of the set from device to host
void opp_download_particle_set(opp_set particles_set, bool force_download) {}

//*******************************************************************************
// Copy all dats of the set from host to device
void opp_upload_particle_set(opp_set particles_set, bool realloc) {}

//*******************************************************************************
void opp_colour_cartesian_mesh(const int ndim, std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors, const int cell_ghosts) {}

//*******************************************************************************
void* opp_host_malloc(size_t size)
{
    return malloc(size);
}

//*******************************************************************************
void* opp_host_realloc(void* ptr, size_t new_size)
{
    return realloc(ptr, new_size);
}

//*******************************************************************************
void opp_host_free(void* ptr)
{
    free(ptr);
}