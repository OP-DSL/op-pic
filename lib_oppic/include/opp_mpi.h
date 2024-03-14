
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

#include <opp_mpi_core.h>

/*******************************************************************************/
void opp_halo_create();
void opp_halo_destroy();

int opp_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args);
void opp_mpi_halo_exchange(oppic_arg *arg, int exec_flag);
void opp_mpi_halo_wait_all(int nargs, oppic_arg *args);
/*******************************************************************************/

void opp_exchange_double_indirect_reductions(oppic_dat dat, opp_reduc_comm reduc_comm);
void opp_complete_double_indirect_reductions(oppic_dat dat);

void cleanSendRecvBuffers(oppic_set set);


inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
    opp_mpi_reduce_double(args, data);
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
    opp_mpi_reduce_int(args, data);
}