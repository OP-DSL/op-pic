
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

#include <oppic_lib_core.h>

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

oppic_set opp_decl_mesh_set(int size, char const *name);

oppic_map opp_decl_mesh_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name);
oppic_dat opp_decl_mesh_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name);

oppic_map oppic_decl_map_txt(oppic_set from, oppic_set to, int dim, const char* file_name, char const *name);
oppic_dat oppic_decl_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name);

oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping = OPP_Map_Default);
oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping = OPP_Map_Default);
oppic_arg opp_get_arg(oppic_dat dat, oppic_access acc, opp_mapping mapping = OPP_Map_Default);
oppic_arg opp_get_arg(oppic_map map, oppic_access acc, opp_mapping mapping = OPP_Map_Default);
oppic_arg opp_get_arg(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping = OPP_Map_Default);

// template <class T> oppic_arg opp_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg opp_get_arg_gbl(double *data, int dim, char const *typ, oppic_access acc);
oppic_arg opp_get_arg_gbl(int *data, int dim, char const *typ, oppic_access acc);
oppic_arg opp_get_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc);

oppic_set opp_decl_part_set(int size, char const *name, oppic_set cells_set);
oppic_set opp_decl_part_set(char const *name, oppic_set cells_set);

oppic_dat opp_decl_part_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index = false);

oppic_dat oppic_decl_particle_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index = false);

template <typename T> 
void opp_decl_const(int dim, T* data, const char* name)
{
   opp_decl_const_impl(dim, sizeof(T), (char*)data, name);
}

void oppic_increase_particle_count(oppic_set particles_set, int num_particles_to_insert);

void opp_inc_part_count_with_distribution(oppic_set particles_set, int num_particles_to_insert, oppic_dat part_dist, bool calc_new = true);

void oppic_reset_num_particles_to_insert(oppic_set set);

void opp_init_particle_move(oppic_set set);
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args);

void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status);

bool opp_finalize_particle_move(oppic_set set);

void oppic_mark_particle_to_remove(oppic_set set, int particle_index);

void oppic_remove_marked_particles_from_set(oppic_set set);
void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove);

void oppic_particle_sort(oppic_set set);

void opp_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix);
void opp_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix);

void oppic_dump_dat(oppic_dat data);

opp_dat opp_fetch_data(opp_dat dat);

void opp_reset_dat(oppic_dat dat, char* val, opp_reset reset = OPP_Reset_All);

opp_move_var opp_get_move_var(int thread = 0);

bool opp_part_check_status(opp_move_var& m, int map0idx, oppic_set set, int particle_index, int& remove_count, int thread = 0);

bool is_double_indirect_reduction(oppic_arg& arg);
void opp_init_double_indirect_reductions(int nargs, oppic_arg *args);
void opp_exchange_double_indirect_reductions(int nargs, oppic_arg *args) ;
void opp_complete_double_indirect_reductions(int nargs, oppic_arg *args);

int opp_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device);
int opp_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args);
void opp_mpi_halo_exchange(oppic_arg *arg, int exec_flag);
void opp_mpi_halo_wait_all(int nargs, oppic_arg *args);

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
void opp_colour_cartesian_mesh(const int ndim, const std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors);
