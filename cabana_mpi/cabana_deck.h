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

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "opp_params.h"

inline OPP_REAL courant_length(OPP_REAL lx, OPP_REAL ly, OPP_REAL lz, size_t nx, size_t ny, size_t nz ) {
    OPP_REAL w0, w1 = 0;
    if( nx>1 ) w0 = nx / lx, w1 += w0*w0;
    if( ny>1 ) w0 = ny / ly, w1 += w0*w0;
    if( nz>1 ) w0 = nz / lz, w1 += w0*w0;
    return sqrt(1 / w1);
}

class Deck {
    public:
    Deck(opp::Params& param) {
        num_steps = param.get<OPP_INT>("num_steps");
        v0 = param.get<OPP_REAL>("init_vel");
        nppc = param.get<OPP_INT>("num_part_per_cell");
        nx = param.get<OPP_INT>("nx");
        ny = param.get<OPP_INT>("ny");
        nz = param.get<OPP_INT>("nz");
        c = param.get<OPP_REAL>("c");
        qsp = param.get<OPP_REAL>("qsp");
        me = param.get<OPP_REAL>("me");

        const OPP_REAL gam = 1.0 / sqrt(1.0 - v0 * v0);
        const OPP_REAL default_grid_len = param.get<OPP_REAL>("default_grid_len");
        len_x_global = default_grid_len;
        len_y_global = 0.628318530717959 * (gam * sqrt(gam));
        len_z_global = default_grid_len;      
        dx = len_x_global / nx;
        dy = len_y_global / ny;
        dz = len_z_global / nz;
        
        n0 = 2.0;
        dt = 0.99 * courant_length(len_x_global, len_y_global, len_z_global, nx, ny, nz) / c;
        
        Npe = n0 * len_x_global * len_y_global * len_z_global;
        Ne = nx * ny * nz * nppc;
        we = Npe / Ne;
        eps = param.get<OPP_REAL>("eps");

        cdt_dx = (c * dt / dx);
        cdt_dy = (c * dt / dy);
        cdt_dz = (c * dt / dz);
        qdt_2mc = qsp * dt / (2 * me * c);
        const OPP_REAL frac = 1.0;
        px = (nx>0) ? (frac * c * dt / dx) : 0;
        py = (ny>0) ? (frac * c * dt / dy) : 0;
        pz = (nz>0) ? (frac * c * dt / dz) : 0;
        acc_coefx = 0.25 / (dy * dz * dt);
        acc_coefy = 0.25 / (dz * dx * dt);
        acc_coefz = 0.25 / (dx * dy * dt);
        dt_eps0 = dt / eps;
    }

    inline void print() {
        opp_printf("Deck", "num_steps: %d ", num_steps);
        opp_printf("Deck", "v0: %2.25lE ", v0);
        opp_printf("Deck", "nppc: %d ", nppc);
        opp_printf("Deck", "nx: %d %d %d", nx, ny, nz);
        opp_printf("Deck", "len_x_global: %2.25lE %2.25lE %2.25lE", len_x_global, len_y_global, len_z_global);
        opp_printf("Deck", "dx: %2.25lE %2.25lE %2.25lE", dx, dy, dz);
        opp_printf("Deck", "n0: %2.25lE ", n0);
        opp_printf("Deck", "dt: %2.25lE ", dt);
        opp_printf("Deck", "c: %2.25lE ", c);
        opp_printf("Deck", "Npe: %2.25lE ", Npe);
        opp_printf("Deck", "Ne: %"PRId64"", Ne);
        opp_printf("Deck", "we: %2.25lE ", we);
        opp_printf("Deck", "eps: %2.25lE ", eps);
        opp_printf("Deck", "qsp: %2.25lE ", qsp);
        opp_printf("Deck", "me: %2.25lE ", me);

        opp_printf("Deck", "cdt_dx: %2.25lE %2.25lE %2.25lE ", cdt_dx, cdt_dy, cdt_dz);
        opp_printf("Deck", "px: %2.25lE %2.25lE %2.25lE ", px, py, pz);
        opp_printf("Deck", "acc_coefx: %2.25lE %2.25lE %2.25lE ", acc_coefx, acc_coefy, acc_coefz);
        opp_printf("Deck", "qdt_2mc: %2.25lE ", qdt_2mc);
        opp_printf("Deck", "dt_eps0: %2.25lE ", dt_eps0);
    }

    OPP_INT num_steps = 0;
    OPP_REAL v0 = 0.0;
    OPP_INT nppc = 0;
    OPP_INT nx = 0;
    OPP_INT ny = 0;
    OPP_INT nz = 0;
    OPP_REAL len_x_global = 0;
    OPP_REAL len_y_global = 0;
    OPP_REAL len_z_global = 0;
    OPP_REAL dx = 0;
    OPP_REAL dy = 0;
    OPP_REAL dz = 0;
    OPP_REAL n0 = 0;
    OPP_REAL dt = 0;
    OPP_REAL c = 0;
    OPP_REAL Npe = 0;
    int64_t Ne = 0;
    OPP_REAL we = 0;
    OPP_REAL eps = 0;
    OPP_REAL qsp = 0;
    OPP_REAL me = 0;

    OPP_REAL cdt_dx = 0;
    OPP_REAL cdt_dy = 0;
    OPP_REAL cdt_dz = 0;
    OPP_REAL qdt_2mc = 0;
    OPP_REAL px = 0;
    OPP_REAL py = 0;
    OPP_REAL pz = 0;
    OPP_REAL acc_coefx = 0;
    OPP_REAL acc_coefy = 0;
    OPP_REAL acc_coefz = 0;
    OPP_REAL dt_eps0 = 0;
};