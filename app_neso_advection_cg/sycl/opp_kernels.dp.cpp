
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

//*********************************************
// AUTO GENERATED CODE
//*********************************************

#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include "opp_sycl.h"

static dpct::constant_memory<OPP_REAL, 1> CONST_extents_d(2);
static dpct::constant_memory<OPP_REAL, 1> CONST_dt_d(1);
static dpct::constant_memory<OPP_REAL, 1> CONST_cell_width_d(1);
static dpct::constant_memory<OPP_INT, 1> CONST_ndimcells_d(2);

static dpct::constant_memory<int, 0> OPP_cells_set_size_d;
int OPP_cells_set_size;

static dpct::constant_memory<int, 0> OPP_comm_iteration_d;

void opp_decl_const_impl(int dim, int size, char *data, const char *name) {
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.in_order_queue();

    if (!strcmp(name, "CONST_extents")) {
        cutilSafeCall(DPCT_CHECK_ERROR(
            q_ct1.memcpy(CONST_extents_d.get_ptr(), data, dim * size).wait()));
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        cutilSafeCall(DPCT_CHECK_ERROR(
            q_ct1.memcpy(CONST_dt_d.get_ptr(), data, dim * size).wait()));
        return;
    }
    if (!strcmp(name, "CONST_cell_width")) {
        cutilSafeCall(DPCT_CHECK_ERROR(
            q_ct1.memcpy(CONST_cell_width_d.get_ptr(), data, dim * size)
                .wait()));
        return;
    }
    if (!strcmp(name, "CONST_ndimcells")) {
        cutilSafeCall(DPCT_CHECK_ERROR(
            q_ct1.memcpy(CONST_ndimcells_d.get_ptr(), data, dim * size)
                .wait()));
        return;
    }

    opp_printf("Error: unknown const name %s", name);
    opp_abort("Error: unknown const name");
}

#include "update_pos_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "verify_kernel_loop.hpp"

