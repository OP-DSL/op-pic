
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

void oppic_init(int argc, char **argv, int diags);
void oppic_exit();

oppic_set oppic_decl_set(int size, char const *name);

oppic_map oppic_decl_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name);

oppic_dat oppic_decl_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name);

oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat(oppic_dat dat, oppic_access acc, bool map_with_cell_index = false);
oppic_arg oppic_arg_dat(oppic_map map, oppic_access acc, bool map_with_cell_index = false);

// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(double *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc);

oppic_set oppic_decl_particle_set(int size, char const *name, oppic_set cells_set);
oppic_set oppic_decl_particle_set(char const *name, oppic_set cells_set);

oppic_dat oppic_decl_particle_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index = false);

void oppic_increase_particle_count(oppic_set particles_set, const int num_particles_to_insert);

void oppic_reset_num_particles_to_insert(oppic_set set);

void oppic_mark_particle_to_remove(oppic_set set, int particle_index);

void oppic_remove_marked_particles_from_set(oppic_set set);
void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove);

void oppic_particle_sort(oppic_set set);

