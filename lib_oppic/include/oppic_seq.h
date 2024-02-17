
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

#include <oppic_lib.h>


template <typename... T, typename... OPARG>
void oppic_par_loop(void (*kernel)(T *...), char const *name, oppic_set set, oppic_iterate_type iter_type,
                 OPARG... arguments) {
    printf("oppic_par_loop %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}

template <typename... T, typename... OPARG>
void oppic_par_loop_particle(void (*kernel)(T *...), char const *name, oppic_set set, oppic_iterate_type iter_type,
                 OPARG... arguments) {
    printf("oppic_par_looppic_particle %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}


inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
    (void)args;
    (void)data;
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
    (void)args;
    (void)data;
}