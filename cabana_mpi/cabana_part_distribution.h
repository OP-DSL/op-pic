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
// USER WRITTEN CODE
//*********************************************

// include only to cabana_misc.cpp
#pragma once

//*************************************************************************************************
inline void enrich_particles_random(const Deck& deck, OPP_INT cell_count, OPP_REAL* pos, OPP_REAL* vel, 
    OPP_REAL* str_mid, OPP_INT* cid, OPP_INT* idx, OPP_REAL* weight, OPP_REAL* cell_pos_ll, OPP_INT rank_part_start) {
    
    printf("Setup[0][0] - enrich_particles_random\n");

    const OPP_INT npart_per_cell = deck.nppc;
    const OPP_REAL const_weight  = deck.we;
    const OPP_REAL init_vel      = deck.v0;
    const OPP_REAL extents[DIM]  = { deck.dx, deck.dy, deck.dz };

    std::mt19937 rng(52234234);
    std::uniform_real_distribution<OPP_REAL> uniform_rng(0.0, 1.0);

    std::vector<OPP_REAL> uniform_dist_vec(npart_per_cell * DIM);
    for (int px = 0; px < npart_per_cell; px++) {
        for (int dimx = 0; dimx < DIM; dimx++) {
            uniform_dist_vec[px * DIM + dimx] = extents[dimx] * uniform_rng(rng);
        }
    }

    // std::string log = "";
    // for (const OPP_REAL& num : uniform_dist_vec)
    //     log += str(num, "%.15f ");
    // opp_printf("W", "%s", log.c_str());

    // Populate the host space with particle data.
    int p_idx = 0;
    for (int cx = 0; cx < cell_count; cx++) {
        for (int px = 0; px < npart_per_cell; px++) {
            
            pos[p_idx * DIM + Dim::x] = uniform_dist_vec[px * DIM + Dim::x] + 
                                                            cell_pos_ll[cx * DIM + Dim::x];
            pos[p_idx * DIM + Dim::y] = uniform_dist_vec[px * DIM + Dim::y] + 
                                                            cell_pos_ll[cx * DIM + Dim::y];
            pos[p_idx * DIM + Dim::z] = uniform_dist_vec[px * DIM + Dim::z] + 
                                                            cell_pos_ll[cx * DIM + Dim::z];

            // pos[p_idx * DIM + Dim::x] = uniform_dist_vec[px * DIM + Dim::x];
            // pos[p_idx * DIM + Dim::y] = uniform_dist_vec[px * DIM + Dim::y];
            // pos[p_idx * DIM + Dim::z] = uniform_dist_vec[px * DIM + Dim::z];

            const int sign = (px % 2) ? 1 : -1; // For two stream
            // const int sign = 1; // For one stream

            vel[p_idx * DIM + Dim::x] = 0.0;
            vel[p_idx * DIM + Dim::y] = sign * init_vel;
            vel[p_idx * DIM + Dim::z] = 0.0;

            str_mid[p_idx * DIM + Dim::x] = 0.0;
            str_mid[p_idx * DIM + Dim::y] = 0.0;
            str_mid[p_idx * DIM + Dim::z] = 0.0;

            cid[p_idx]    = cx;
            idx[p_idx]    = (cx * npart_per_cell + px + rank_part_start); // this might not exactly match with MPI versions
            weight[p_idx] = const_weight;

            p_idx++;
        }
    }
}

//*************************************************************************************************
inline void enrich_particles_two_stream(const Deck& deck, OPP_INT cell_count, OPP_REAL* pos,  
                    OPP_REAL* vel, OPP_REAL* str_mid, OPP_INT* cid, OPP_INT* idx, OPP_REAL* weight,
                    OPP_REAL* cell_pos_ll, OPP_INT rank_part_start, OPP_INT* global_cids) {
    
    printf("Setup[0][0] - enrich_particles_two_stream\n");

    const OPP_INT npart_per_cell = deck.nppc;
    const OPP_REAL const_weight  = deck.we;
    const OPP_REAL v0            = deck.v0;
    const OPP_INT ny             = deck.ny;
    const OPP_REAL dxp           = (2.0 / npart_per_cell);
    const OPP_REAL n_particles   = (npart_per_cell * cell_count);

    // Populate the host space with particle data.
    for (int p_idx = 0; p_idx < n_particles; p_idx++) {

        int sign =  -1;
        const size_t pi2 = (rank_part_start + p_idx);
        const size_t pi = ((pi2) / 2);
        if (pi2 % 2 == 0)
            sign = 1;

        const int local_cid = (2 * pi / npart_per_cell);
        const int global_cid = global_cids[local_cid];

        const int pic = (2 * pi) % npart_per_cell; //Every 2 particles have the same "pic".
        const double x = pic * dxp + 0.5 * dxp - 1.0;

        // Initialize velocity.(each cell length is 2)
        const double gam = 1.0 / sqrt(1.0 - v0 * v0);

        const double na = 0.0001 * sin(2.0 * 3.1415926 * ((x + 1.0 + global_cid * 2) / (2 * ny)));
        
        pos[p_idx * DIM + Dim::x] = 0.0;
        pos[p_idx * DIM + Dim::y] = x;
        pos[p_idx * DIM + Dim::z] = 0.0;

        vel[p_idx * DIM + Dim::x] = sign * v0 * gam * (1.0 + na * sign);
        vel[p_idx * DIM + Dim::y] = 0;
        vel[p_idx * DIM + Dim::z] = 0;

        str_mid[p_idx * DIM + Dim::x] = 0.0;
        str_mid[p_idx * DIM + Dim::y] = 0.0;
        str_mid[p_idx * DIM + Dim::z] = 0.0;

        cid[p_idx]    = local_cid;
        idx[p_idx]    = (p_idx + rank_part_start); // this might not exactly match with MPI versions
        weight[p_idx] = const_weight;  
    }
}