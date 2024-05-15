
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

#pragma once

#include <opp_lib_core.h>

const double opp_zero_double16[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
const int opp_zero_int16[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

extern std::vector<int> part_remove_count_per_thr;

struct opp_move_var
{
    bool inside_cell = true;
    bool iteration_one = true;
    opp_move_status move_status = OPP_MOVE_DONE;
};

void opp_init(int argc, char **argv);
void opp_exit();

opp_set opp_decl_set(int size, char const *name);

opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name);

opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index = false);

opp_dat opp_decl_mesh_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name);

opp_map opp_decl_map_txt(opp_set from, opp_set to, int dim, const char* file_name, char const *name);
opp_dat opp_decl_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name);

opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_dat p2c_map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_dat p2c_map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_dat dat, opp_dat p2c_map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_dat dat, opp_access acc, bool offset = true);

opp_arg opp_arg_dat(opp_map data_map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_map data_map, opp_dat p2c_map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_access acc, bool offset = true);
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_dat p2c_map, opp_access acc, bool offset = true);

// template <class T> opp_arg opp_arg_gbl(T *data, int dim, char const *typ, opp_access acc);
opp_arg opp_arg_gbl(double *data, int dim, char const *typ, opp_access acc);
opp_arg opp_arg_gbl(int *data, int dim, char const *typ, opp_access acc);
opp_arg opp_arg_gbl(const bool *data, int dim, char const *typ, opp_access acc);

opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set);
opp_set opp_decl_particle_set(char const *name, opp_set cells_set);

opp_dat opp_decl_part_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index = false);

opp_dat opp_decl_particle_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index = false);

template <typename T> 
void opp_decl_const(int dim, T* data, const char* name)
{
    opp_decl_const_impl(dim, sizeof(T), (char*)data, name);
}

void opp_increase_particle_count(opp_set particles_set, int num_particles_to_insert);

void opp_inc_part_count_with_distribution(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist, bool calc_new = true);

void opp_reset_num_particles_to_insert(opp_set set);

void opp_init_particle_move(opp_set set);
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args);

void opp_mark_particle_to_move(opp_set set, int particle_index, int move_status);

bool opp_finalize_particle_move(opp_set set);

void opp_mark_particle_to_remove(opp_set set, int particle_index);

void opp_remove_marked_particles_from_set(opp_set set);
void opp_remove_marked_particles_from_set(opp_set set, std::vector<int>& idx_to_remove);

void opp_particle_sort(opp_set set);

void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix);
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix);

void opp_dump_dat(opp_dat data);

opp_dat opp_fetch_data(opp_dat dat);

void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset = OPP_Reset_All);

template <typename T> 
void opp_reset_dat(opp_dat dat, T* val, opp_reset reset = OPP_Reset_All)
{
    opp_reset_dat_impl(dat, (char*)val, reset);
}

opp_move_var opp_get_move_var(int thread = 0);

bool opp_part_check_status(opp_move_var& m, int map0idx, opp_set set, int particle_index, int& remove_count, int thread = 0);

bool is_double_indirect_reduction(opp_arg& arg);
void opp_init_double_indirect_reductions(int nargs, opp_arg *args);
void opp_exchange_double_indirect_reductions(int nargs, opp_arg *args) ;
void opp_complete_double_indirect_reductions(int nargs, opp_arg *args);

int opp_mpi_halo_exchanges_grouped(opp_set set, int nargs, opp_arg *args, DeviceType device);
int opp_mpi_halo_exchanges(opp_set set, int nargs, opp_arg *args);
void opp_mpi_halo_exchange(opp_arg *arg, int exec_flag);
void opp_mpi_halo_wait_all(int nargs, opp_arg *args);

// Copy a dat from host to device
void opp_upload_dat(opp_dat dat);

// Copy a dat from device to host
void opp_download_dat(opp_dat dat);

// Copy all dats of the set from device to host
void opp_download_particle_set(opp_set particles_set, bool force_download = false);

// Copy all dats of the set from host to device
void opp_upload_particle_set(opp_set particles_set, bool realloc = false);

//*************************************************************************************************
// ndim : can be 2 or 3
// cell_counts : cell_counts in each direction
// cell_index : cell_index dat which holds global numbering
// cell_colors : local cell_colors dat to colour with most appropriate MPI rank
void opp_colour_cartesian_mesh(const int ndim, std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors, const int cell_ghosts = 0);

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif