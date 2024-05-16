
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

#include "opp_seq.h"

template <typename... T, typename... OPARG>
void opp_par_loop(void (*kernel)(T *...), char const *name, opp_set set, opp_iterate_type iter_type,
                 OPARG... arguments) {
    opp_printf("opp_par_loop", "kernel %s iterate %s", name, (iter_type == OPP_ITERATE_ALL) ? "ALL" : "INJECTED");
}

template <typename... T, typename... OPARG>
void opp_particle_move(void (*kernel)(T *...), char const *name, opp_set set, opp_map c2c_map, opp_dat p2c_map,
                 OPARG... arguments) {
    opp_printf("opp_particle_move", "kernel %s", name);
}

inline void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    opp_printf("opp_decl_const", "name %s dim %d size %d", name, dim, size);
}

inline void opp_init_direct_hop(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& bounding_box)
{
    opp_printf("opp_init_direct_hop", "grid_spacing %lf dim %d c_gbl_id %s", grid_spacing, dim, c_gbl_id->name);
}