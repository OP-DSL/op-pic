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

#include "cabana_defs.h"

using namespace opp;

#include "cabana_misc.cpp"

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

        OPP_INT max_iter = opp_params->get<OPP_INT>("max_iter");   
        OPP_REAL dt      = opp_params->get<OPP_REAL>("dt");
        OPP_REAL qsp     = opp_params->get<OPP_REAL>("qsp");
        OPP_REAL qdt_2mc = cabana_get_qdt_2mc();
        OPP_REAL dt_eps0 = cabana_get_dt_eps0();

        std::shared_ptr<DataPointers> m = LoadData();

        opp_set cell_set        = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_map cell_cell_map   = opp_decl_mesh_map(cell_set, cell_set, NEIGHBOURS, m->cell_cell_map, "c_c_map");
        opp_dat cell_index      = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_index,  "c_index");
        opp_dat cell_e          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_e,      "c_e");  
        opp_dat cell_b          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_b,      "c_b");        
        opp_dat cell_j          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_j,      "c_j");
        opp_dat cell_acc        = opp_decl_mesh_dat(cell_set, ACC_LEN,    DT_REAL, m->c_acc,    "c_acc"); 
        opp_dat cell_interp     = opp_decl_mesh_dat(cell_set, INTERP_LEN, DT_REAL, m->c_interp, "c_interp"); 
        opp_dat cell_pos_ll     = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_pos_ll, "c_pos_ll");       
        opp_dat cell_colors     = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_index,  "c_colors"); // init dummy

        opp_set part_set        = opp_decl_part_set("particles", cell_set);  // Zero particles, inject after partitioning
        opp_dat part_index      = opp_decl_part_dat(part_set, ONE, DT_INT,  nullptr, "p_index"); // Unused
        opp_dat part_pos        = opp_decl_part_dat(part_set, DIM, DT_REAL, nullptr, "p_pos");
        opp_dat part_vel        = opp_decl_part_dat(part_set, DIM, DT_REAL, nullptr, "p_vel");    
        opp_dat part_streak_mid = opp_decl_part_dat(part_set, DIM, DT_REAL, nullptr, "p_streak_mid");
        opp_dat part_weight     = opp_decl_part_dat(part_set, ONE, DT_REAL, nullptr, "p_weight");
        opp_dat part_mesh_rel   = opp_decl_part_dat(part_set, ONE, DT_INT,  nullptr, "p_mesh_rel", true);

        opp_decl_const<OPP_INT>(DIM,  cabana_get_cells_per_dim().data(), "CONST_c_per_dim");
        opp_decl_const<OPP_REAL>(ONE, &dt,                               "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &qsp,                              "CONST_qsp");
        opp_decl_const<OPP_REAL>(DIM, cabana_get_cdt_d().data(),         "CONST_cdt_d");
        opp_decl_const<OPP_REAL>(DIM, cabana_get_p().data(),             "CONST_p");
        opp_decl_const<OPP_REAL>(ONE, &qdt_2mc,                          "CONST_qdt_2mc");
        opp_decl_const<OPP_REAL>(ONE, &dt_eps0,                          "CONST_dt_eps0");
        opp_decl_const<OPP_REAL>(DIM, cabana_get_acc_coef().data(),      "CONST_acc_coef");

        m->DeleteValues();

        // ideally opp_colour_cartesian_mesh is not required for non-mpi runs
        opp_colour_cartesian_mesh(DIM, cabana_get_cells_per_dim(), cell_index, cell_colors);

#ifdef USE_MPI
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif

        init_particles(part_index, part_pos, part_vel, part_streak_mid, part_weight, part_mesh_rel, 
                        cell_pos_ll);

        opp_printf("Setup", "Cells[%d] Particles[%d]", cell_set->size, part_set->size);

        opp_profiler->end("Setup");

        opp_profiler->start("MainLoop");

        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            if (OP_DEBUG && OPP_rank == OPP_ROOT) 
                opp_printf("Main", "Starting main loop iteration %d ****", OPP_main_loop_iter);

            opp_loop_all__interpolate_mesh_fields(
                cell_set,
                opp_get_arg(cell_e,                                  OP_READ),
                opp_get_arg(cell_b,                                  OP_READ),
                opp_get_arg(cell_e, CellMap::xu_y_z,  cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::x_yu_z,  cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::x_y_zu,  cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::xu_yu_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::x_yu_zu, cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::xu_y_zu, cell_cell_map, OP_READ),
                opp_get_arg(cell_b, CellMap::xu_y_z,  cell_cell_map, OP_READ),
                opp_get_arg(cell_b, CellMap::x_yu_z,  cell_cell_map, OP_READ),
                opp_get_arg(cell_b, CellMap::x_y_zu,  cell_cell_map, OP_READ),
                opp_get_arg(cell_interp,                             OP_WRITE)
            );

            opp_reset_dat(
                cell_acc, 
                (char*)opp_zero_double16);

            opp_particle_mover__Move(
                part_set,
                opp_get_arg(part_mesh_rel,   OP_RW),
                opp_get_arg(part_vel,        OP_RW),
                opp_get_arg(part_pos,        OP_RW),
                opp_get_arg(part_streak_mid, OP_RW),
                opp_get_arg(part_weight,     OP_READ),
                opp_get_arg(cell_interp,     OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_acc,        OP_INC,  OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_cell_map,   OP_READ, OPP_Map_from_Mesh_Rel)
            );

            opp_loop_all__accumulate_current_to_cells(
                cell_set,
                opp_get_arg(cell_j,                                    OP_WRITE),
                opp_get_arg(cell_acc,                                  OP_READ),
                opp_get_arg(cell_acc, CellMap::xd_y_z , cell_cell_map, OP_READ),
                opp_get_arg(cell_acc, CellMap::x_yd_z , cell_cell_map, OP_READ),
                opp_get_arg(cell_acc, CellMap::x_y_zd , cell_cell_map, OP_READ),
                opp_get_arg(cell_acc, CellMap::xd_yd_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_acc, CellMap::x_yd_zd, cell_cell_map, OP_READ),
                opp_get_arg(cell_acc, CellMap::xd_y_zd, cell_cell_map, OP_READ) 
            );

            // Leap frog method
            opp_loop_all__half_advance_b(
                cell_set,
                opp_get_arg(cell_e, CellMap::xu_y_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::x_yu_z, cell_cell_map, OP_READ), 
                opp_get_arg(cell_e, CellMap::x_y_zu, cell_cell_map, OP_READ), 
                opp_get_arg(cell_e,                                 OP_READ),
                opp_get_arg(cell_b,                                 OP_INC)
            );

            opp_loop_all__advance_e(
                cell_set,
                opp_get_arg(cell_b, CellMap::xd_y_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_b, CellMap::x_yd_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_b, CellMap::x_y_zd, cell_cell_map, OP_READ),
                opp_get_arg(cell_b,                                 OP_READ),
                opp_get_arg(cell_j,                                 OP_READ),
                opp_get_arg(cell_e,                                 OP_INC) 
            );

            opp_loop_all__half_advance_b(
                cell_set,
                opp_get_arg(cell_e, CellMap::xu_y_z, cell_cell_map, OP_READ),
                opp_get_arg(cell_e, CellMap::x_yu_z, cell_cell_map, OP_READ), 
                opp_get_arg(cell_e, CellMap::x_y_zu, cell_cell_map, OP_READ), 
                opp_get_arg(cell_e,                                 OP_READ),
                opp_get_arg(cell_b,                                 OP_INC) 
            );

            std::string log = ""; // TODO : print some unseful information to verify
            if (opp_params->get<OPP_BOOL>("print_final"))
            {
                OPP_REAL max_j = 0.0, max_e = 0.0, max_b = 0.0;

                opp_loop_all__GetFinalMaxValues( // plan is to get only the x values reduced here
                    cell_set,
                    opp_get_arg(cell_j, OP_READ),
                    opp_get_arg_gbl(&max_j, 1, "double", OP_MAX),
                    opp_get_arg(cell_e, OP_READ),
                    opp_get_arg_gbl(&max_e, 1, "double", OP_MAX),
                    opp_get_arg(cell_b, OP_READ),
                    opp_get_arg_gbl(&max_b, 1, "double", OP_MAX)
                );

                log += str(max_j, " max_j: %2.15lE");
                log += str(max_e, " max_e: %2.15lE");
                log += str(max_b, " max_b: %2.15lE");
            }
            if (OPP_rank == OPP_ROOT) 
                opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
        opp_profiler->end("MainLoop");
        
        if (OPP_rank == OPP_ROOT) 
            opp_printf("Main", "Main loop completed after %d iterations ****", max_iter);
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(cell_v_nodes_map  , f.c_str(), "cell_v_nodes_map.dat");
// opp_print_dat_to_txtfile(node_charge_den, f.c_str(), "node_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(cell_shape_deriv, "cell_shape_deriv.dat");


//             if (false)
//             { 
//                 std::string midString = "A_ACC_";  
// #ifdef USE_MPI  
//                 std::string f1 = midString + std::to_string(OPP_main_loop_iter) + "_cell_e.dat";
//                 std::string f2 = midString + std::to_string(OPP_main_loop_iter) + "_cell_b.dat";
//                 std::string f3 = midString + std::to_string(OPP_main_loop_iter) + "_cell_j.dat";
//                 std::string f4 = midString + std::to_string(OPP_main_loop_iter) + "_cell_interp.dat";
//                 std::string f5 = midString + std::to_string(OPP_main_loop_iter) + "_cell_acc.dat";

//                 opp_mpi_print_dat_to_txtfile(cell_e, f1.c_str());
//                 opp_mpi_print_dat_to_txtfile(cell_b, f2.c_str());
//                 opp_mpi_print_dat_to_txtfile(cell_j, f3.c_str());
//                 opp_mpi_print_dat_to_txtfile(cell_interp, f4.c_str());
//                 opp_mpi_print_dat_to_txtfile(cell_acc, f5.c_str());
// #else           
//                 std::string f = std::string("F_") + midString + std::to_string(OPP_main_loop_iter);
                
//                 opp_print_dat_to_txtfile(cell_e, f.c_str(), "cell_e.dat");
//                 opp_print_dat_to_txtfile(cell_b, f.c_str(), "cell_b.dat");
//                 opp_print_dat_to_txtfile(cell_j, f.c_str(), "cell_j.dat");
//                 opp_print_dat_to_txtfile(cell_interp, f.c_str(), "cell_interp.dat");
//                 opp_print_dat_to_txtfile(cell_acc, f.c_str(), "cell_acc.dat");
// #endif
//             }