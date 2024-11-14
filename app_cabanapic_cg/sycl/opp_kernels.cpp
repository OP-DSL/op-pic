
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

#include "opp_sycl.h"

OPP_REAL* CONST_dt_s = nullptr;
OPP_REAL* CONST_qsp_s = nullptr;
OPP_REAL* CONST_cdt_d_s = nullptr;
OPP_REAL* CONST_p_s = nullptr;
OPP_REAL* CONST_qdt_2mc_s = nullptr;
OPP_REAL* CONST_dt_eps0_s = nullptr;
OPP_REAL* CONST_acc_coef_s = nullptr;
    
OPP_INT* cells_set_size_s = nullptr;
OPP_INT cells_set_size = -1;

OPP_INT* comm_iteration_s = nullptr;
OPP_INT comm_iteration = -1;

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_dt")) {
        opp_register_const<OPP_REAL>(CONST_dt_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_dt_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_qsp")) {
        opp_register_const<OPP_REAL>(CONST_qsp_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_qsp_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_cdt_d")) {
        opp_register_const<OPP_REAL>(CONST_cdt_d_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_cdt_d_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_p")) {
        opp_register_const<OPP_REAL>(CONST_p_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_p_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_qdt_2mc")) {
        opp_register_const<OPP_REAL>(CONST_qdt_2mc_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_qdt_2mc_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_dt_eps0")) {
        opp_register_const<OPP_REAL>(CONST_dt_eps0_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_dt_eps0_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_acc_coef")) {
        opp_register_const<OPP_REAL>(CONST_acc_coef_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_acc_coef_s, (OPP_REAL*)data, dim);
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "interpolate_mesh_fields_kernel_loop.hpp"

#include "move_deposit_kernel_loop.hpp"

#include "accumulate_current_to_cells_kernel_loop.hpp"

#include "half_advance_b_kernel_loop.hpp"

#include "update_ghosts_B_kernel_loop.hpp"

#include "update_ghosts_kernel_loop.hpp"

#include "advance_e_kernel_loop.hpp"

#include "compute_energy_kernel_loop.hpp"

#include "get_max_x_values_kernel_loop.hpp"

