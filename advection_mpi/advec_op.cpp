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

#include "advec_defs.h"

using namespace opp;

#include "advec_misc.cpp"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    {
        opp_profiler->start("Setup");

        OPP_INT max_iter          = opp_params->get<OPP_INT>("max_iter");   
        OPP_REAL const_dt         = opp_params->get<OPP_REAL>("dt");
        OPP_REAL cell_width       = opp_params->get<OPP_REAL>("cell_width");
        OPP_REAL const_extents[2] = {opp_params->get<OPP_INT>("nx")*cell_width, opp_params->get<OPP_INT>("ny")*cell_width};

        std::shared_ptr<DataPointers> m = LoadData();

        opp_set cell_set      = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_map cell_cell_map = opp_decl_mesh_map(cell_set, cell_set, NEIGHBOURS, m->cell_cell_map, "c_c_map");
        opp_dat cell_index    = opp_decl_mesh_dat(cell_set, ONE, DT_INT,  m->c_index,  "c_index");
        opp_dat cell_pos_ll   = opp_decl_mesh_dat(cell_set, DIM, DT_REAL, m->c_pos_ll, "c_pos_ll");       
        opp_dat cell_colors   = opp_decl_mesh_dat(cell_set, ONE, DT_INT,  m->c_colors, "c_colors"); // used only with MPI

        opp_set part_set      = opp_decl_part_set("particles", cell_set); // Zero particles, inject after partitioning
        opp_dat part_index    = opp_decl_part_dat(part_set, ONE, DT_INT,  nullptr, "p_index");
        opp_dat part_pos      = opp_decl_part_dat(part_set, DIM, DT_REAL, nullptr, "p_pos");
        opp_dat part_vel      = opp_decl_part_dat(part_set, DIM, DT_REAL, nullptr, "p_vel");    
        opp_dat part_mesh_rel = opp_decl_part_dat(part_set, ONE, DT_INT,  nullptr, "p_mesh_rel", true);

        opp_decl_const<OPP_REAL>(TWO, const_extents, "CONST_extents");
        opp_decl_const<OPP_REAL>(ONE, &const_dt,     "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &cell_width,   "CONST_cell_width");

        m->DeleteValues();
        
        // opp_print_map_to_txtfile(cell_cell_map, "INIT_B", "cell_cell_map.dat");
#ifdef USE_MPI
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif
        // opp_print_map_to_txtfile(cell_cell_map, "INIT_A", "cell_cell_map.dat");

        init_particles(part_index, part_pos, part_vel, part_mesh_rel, cell_pos_ll);

        opp_printf("Setup", "Cells[%d] Particles[%d]", cell_set->size, part_set->size);

        opp_profiler->end("Setup");

        opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++) {
        
            opp_particle_mover__Move(
                part_set,
                opp_get_arg(part_mesh_rel, OP_RW),
                opp_get_arg(part_vel,      OP_RW),
                opp_get_arg(part_pos,      OP_RW),
                opp_get_arg(cell_pos_ll,   OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_cell_map, OP_READ, OPP_Map_from_Mesh_Rel)
            );

            if (OPP_rank == OPP_ROOT) opp_printf("Main", "ts: %d ****", OPP_main_loop_iter);
        }
        opp_profiler->end("MainLoop");
        
        if (OPP_rank == OPP_ROOT) 
            opp_printf("Main", "Main loop completed after %d iterations ****", max_iter);

        // opp_print_dat_to_txtfile(part_mesh_rel, "FINAL", "part_mesh_rel.dat");
        // opp_print_dat_to_txtfile(part_pos, "FINAL", "part_pos.dat");
    }

    opp_exit();

    return 0;
}
