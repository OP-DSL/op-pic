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

// make clean seq
// bin/seq config/cabana.param

#include "cabana_defs.h"

using namespace opp;

#include "cabana_misc.cpp"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << 
        "/ext-home/zl/phd/OP-PIC/scripts/configs/coarse_1.param" << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    {
        opp_profiler->start("Setup");

        int max_iter    = opp_params->get<OPP_INT>("max_iter");   
        int const_nx    = opp_params->get<OPP_INT>("nx");
        int const_ny    = opp_params->get<OPP_INT>("ny");
        int const_nz    = opp_params->get<OPP_INT>("nz");
        int const_ng    = opp_params->get<OPP_INT>("ng");
        std::string log = "";

        std::shared_ptr<DataPointers> g_m, m; // g_m - global mesh, m - local mesh
        g_m = std::make_shared<DataPointers>();
        
        if (OPP_rank == OPP_ROOT)
        {
            g_m = LoadData();
        } 
        DistributeDataOverRanks(g_m, m);

        opp_set cell_set        = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_set part_set        = opp_decl_part_set(m->n_particles, "particles", cell_set); 

        // opp_map acc_cell_map    = opp_decl_mesh_map(cell_set, cell_set, 2*DIM, m->acc_cell_map,    "acc_cell_map");   
        // opp_map interp_cell_map = opp_decl_mesh_map(cell_set, cell_set, 2*DIM, m->interp_cell_map, "interp_cell_map"); 
        // opp_map adv_b_cell_map  = opp_decl_mesh_map(cell_set, cell_set, DIM,   m->adv_b_cell_map,  "adv_b_cell_map"); 
        // opp_map adv_e_cell_map  = opp_decl_mesh_map(cell_set, cell_set, DIM,   m->adv_e_cell_map,  "adv_e_cell_map"); 
        // opp_map move_cell_map   = opp_decl_mesh_map(cell_set, cell_set, 2*DIM, m->move_cell_map,   "move_cell_map");  
        
        opp_map cell_cell_map   = opp_decl_mesh_map(cell_set, cell_set, NEIGHBOUR_CELLS, m->cell_cell_map, "cell_cell_map");

        opp_dat cell_index      = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_index,     "c_index");
        opp_dat cell_e          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_e,          "c_e");  
        opp_dat cell_b          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_b,          "c_b");        
        opp_dat cell_j          = opp_decl_mesh_dat(cell_set, DIM,        DT_REAL, m->c_j,          "c_j");
        opp_dat cell_acc        = opp_decl_mesh_dat(cell_set, 12,         DT_REAL, m->c_acc,        "c_acc"); 
        opp_dat cell_interp     = opp_decl_mesh_dat(cell_set, INTERP_LEN, DT_REAL, m->c_interp,     "c_interp"); 
        opp_dat cell_ghost      = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_ghost,      "c_ghost");
        opp_dat cell_iter_adv_e = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_iter_adv_e, "c_adv_e");
        opp_dat cell_iter_acc   = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_iter_acc,   "c_iter_acc");
#ifdef USE_MPI        
        opp_dat cell_colors     = opp_decl_mesh_dat(cell_set, ONE,        DT_INT,  m->c_iter_acc,   "c_colors"); // init dummy
#endif

        opp_dat part_index      = opp_decl_part_dat(part_set, ONE, DT_INT,  m->p_index,  "p_index");
        opp_dat part_pos        = opp_decl_part_dat(part_set, DIM, DT_REAL, m->p_pos,    "p_pos");
        opp_dat part_vel        = opp_decl_part_dat(part_set, DIM, DT_REAL, m->p_vel,    "p_vel");    
        opp_dat part_streak_mid = opp_decl_part_dat(part_set, DIM, DT_REAL, m->p_st_mid, "p_streak_mid");
        opp_dat part_weight     = opp_decl_part_dat(part_set, ONE, DT_REAL, m->p_weight, "p_weight");
        opp_dat part_mesh_rel   = opp_decl_part_dat(part_set, ONE, DT_INT,  m->p_cid,    "p_mesh_rel", true);
        
        opp_decl_const<OPP_INT>(ONE, &const_nx,           "CONST_nx");
        opp_decl_const<OPP_INT>(ONE, &const_ny,           "CONST_ny");
        opp_decl_const<OPP_INT>(ONE, &const_nz,           "CONST_nz");
        opp_decl_const<OPP_INT>(ONE, &const_ng,           "CONST_ng");

        m->DeleteValues();

// opp_print_dat_to_txtfile(cell_index, "INIT", "cell_index_bef_p.dat");
// opp_print_map_to_txtfile(cell_cell_map, "INIT", "cell_cell_map_bef_p.dat");

#ifdef USE_MPI
        genColoursForBlockPartition(cell_colors, cell_ghost, cell_index);

opp_print_dat_to_txtfile(cell_index, "AA", "cell_index.dat");
opp_print_dat_to_txtfile(cell_colors, "AA", "cell_colors.dat");

        // opp_partition(std::string("PARMETIS_KWAY"), cell_set, cell_v_nodes_map);
        // opp_partition(std::string("PARMETIS_GEOM"), iface_set, nullptr, iface_n_pos);
        // opp_partition(std::string("EXTERNAL"), node_set, nullptr, node_colors);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);

opp_printf("Main", "Cells set size after Partitioning %d ****", cell_set->size);

        // updatePartsWithLocalCellIndices(cell_index, part_mesh_rel, part_index);

opp_print_dat_to_txtfile(cell_index, "INIT", "cell_index.dat");
opp_print_dat_to_txtfile(part_mesh_rel, "INIT", "part_mesh_rel.dat");
opp_print_dat_to_txtfile(part_index, "INIT", "part_index.dat");

// opp_print_dat_to_txtfile(cell_index, "INIT", "cell_index_aft_p.dat");
// opp_print_map_to_txtfile(cell_cell_map, "INIT", "cell_cell_map_aft_p.dat");

        MPI_Barrier(MPI_COMM_WORLD);
#endif

opp_printf("Main", "%d number of particles before adding", part_set->size);

init_particles(part_index, part_pos, part_vel, part_streak_mid, 
    part_weight, part_mesh_rel, cell_index, cell_ghost);

opp_printf("Main", "%d number of particles after adding", part_set->size);

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
                opp_get_arg(cell_ghost,                              OP_READ),
                opp_get_arg(cell_e, CellMap::xu_y_z,  cell_cell_map, OP_READ),      // 0, interp_cell_map,    CellMap::interp_0, cell_cell_map
                opp_get_arg(cell_e, CellMap::x_yu_z,  cell_cell_map, OP_READ),      // 1, interp_cell_map,    CellMap::interp_1, cell_cell_map
                opp_get_arg(cell_e, CellMap::x_y_zu,  cell_cell_map, OP_READ),      // 2, interp_cell_map,    CellMap::interp_2, cell_cell_map
                opp_get_arg(cell_e, CellMap::xu_yu_z, cell_cell_map, OP_READ),      // 3, interp_cell_map,    CellMap::interp_3, cell_cell_map
                opp_get_arg(cell_e, CellMap::x_yu_zu, cell_cell_map, OP_READ),      // 4, interp_cell_map,    CellMap::interp_4, cell_cell_map
                opp_get_arg(cell_e, CellMap::xu_y_zu, cell_cell_map, OP_READ),      // 5, interp_cell_map,    CellMap::interp_5, cell_cell_map
                opp_get_arg(cell_b, CellMap::xu_y_z,  cell_cell_map, OP_READ),      // 0, interp_cell_map,    CellMap::interp_0, cell_cell_map
                opp_get_arg(cell_b, CellMap::x_yu_z,  cell_cell_map, OP_READ),      // 1, interp_cell_map,    CellMap::interp_1, cell_cell_map
                opp_get_arg(cell_b, CellMap::x_y_zu,  cell_cell_map, OP_READ),      // 2, interp_cell_map,    CellMap::interp_2, cell_cell_map
                opp_get_arg(cell_interp,                             OP_WRITE)
            );

            if (false)
            { 
                std::string midString = "A_IN_";  
#ifdef USE_MPI  
                std::string f1 = midString + std::to_string(OPP_main_loop_iter) + "_cell_e.dat";
                std::string f2 = midString + std::to_string(OPP_main_loop_iter) + "_cell_b.dat";
                std::string f3 = midString + std::to_string(OPP_main_loop_iter) + "_cell_j.dat";
                std::string f4 = midString + std::to_string(OPP_main_loop_iter) + "_cell_interp.dat";
                std::string f5 = midString + std::to_string(OPP_main_loop_iter) + "_cell_acc.dat";

                opp_mpi_print_dat_to_txtfile(cell_e, f1.c_str());
                opp_mpi_print_dat_to_txtfile(cell_b, f2.c_str());
                opp_mpi_print_dat_to_txtfile(cell_j, f3.c_str());
                opp_mpi_print_dat_to_txtfile(cell_interp, f4.c_str());
                opp_mpi_print_dat_to_txtfile(cell_acc, f5.c_str());
#else           
                std::string f = std::string("F_") + midString + std::to_string(OPP_main_loop_iter);
                
                opp_print_dat_to_txtfile(cell_e, f.c_str(), "cell_e.dat");
                opp_print_dat_to_txtfile(cell_b, f.c_str(), "cell_b.dat");
                opp_print_dat_to_txtfile(cell_j, f.c_str(), "cell_j.dat");
                opp_print_dat_to_txtfile(cell_interp, f.c_str(), "cell_interp.dat");
                opp_print_dat_to_txtfile(cell_acc, f.c_str(), "cell_acc.dat");
#endif
            } // No Change here in iter 0

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
                opp_get_arg(cell_acc,        OP_INC,  OPP_Map_from_Mesh_Rel)
            );

            if (false)
            { 
                std::string midString = "A_MV_";  
#ifdef USE_MPI  
                std::string f1 = midString + std::to_string(OPP_main_loop_iter) + "_cell_e.dat";
                std::string f2 = midString + std::to_string(OPP_main_loop_iter) + "_cell_b.dat";
                std::string f3 = midString + std::to_string(OPP_main_loop_iter) + "_cell_j.dat";
                std::string f4 = midString + std::to_string(OPP_main_loop_iter) + "_cell_interp.dat";
                std::string f5 = midString + std::to_string(OPP_main_loop_iter) + "_cell_acc.dat";

#ifdef USE_MPI
    int nargs = 1;
    oppic_arg args[nargs];
    args[0] = opp_get_arg(cell_acc,    OP_RW);
    args[0].idx = 2; // HACK to forcefully make halos to download

    opp_mpi_halo_exchanges_grouped(cell_acc->set, nargs, args, Device_CPU);
    opp_mpi_halo_wait_all(nargs, args);
#endif

                opp_mpi_print_dat_to_txtfile(cell_e, f1.c_str());
                opp_mpi_print_dat_to_txtfile(cell_b, f2.c_str());
                opp_mpi_print_dat_to_txtfile(cell_j, f3.c_str());
                opp_mpi_print_dat_to_txtfile(cell_interp, f4.c_str());
                opp_mpi_print_dat_to_txtfile(cell_acc, f5.c_str());
#else           
                std::string f = std::string("F_") + midString + std::to_string(OPP_main_loop_iter);
                
                opp_print_dat_to_txtfile(cell_e, f.c_str(), "cell_e.dat");
                opp_print_dat_to_txtfile(cell_b, f.c_str(), "cell_b.dat");
                opp_print_dat_to_txtfile(cell_j, f.c_str(), "cell_j.dat");
                opp_print_dat_to_txtfile(cell_interp, f.c_str(), "cell_interp.dat");
                opp_print_dat_to_txtfile(cell_acc, f.c_str(), "cell_acc.dat");
#endif
            } // cell_acc Change in iter 0

            opp_loop_all__accumulate_current_to_cells(
                cell_set,
                opp_get_arg(cell_iter_acc,                             OP_READ),
                opp_get_arg(cell_j,                                    OP_WRITE),
                opp_get_arg(cell_acc,                                  OP_READ),
                opp_get_arg(cell_acc, CellMap::xd_y_z , cell_cell_map, OP_READ),    // 0, acc_cell_map,
                opp_get_arg(cell_acc, CellMap::x_yd_z , cell_cell_map, OP_READ),    // 1, acc_cell_map,
                opp_get_arg(cell_acc, CellMap::x_y_zd , cell_cell_map, OP_READ),    // 2, acc_cell_map,
                opp_get_arg(cell_acc, CellMap::xd_yd_z, cell_cell_map, OP_READ),    // 3, acc_cell_map,
                opp_get_arg(cell_acc, CellMap::x_yd_zd, cell_cell_map, OP_READ),    // 4, acc_cell_map,
                opp_get_arg(cell_acc, CellMap::xd_y_zd, cell_cell_map, OP_READ)     // 5, acc_cell_map,
            );

            if (false)
            { 
                std::string midString = "A_ACC_";  
#ifdef USE_MPI  
                std::string f1 = midString + std::to_string(OPP_main_loop_iter) + "_cell_e.dat";
                std::string f2 = midString + std::to_string(OPP_main_loop_iter) + "_cell_b.dat";
                std::string f3 = midString + std::to_string(OPP_main_loop_iter) + "_cell_j.dat";
                std::string f4 = midString + std::to_string(OPP_main_loop_iter) + "_cell_interp.dat";
                std::string f5 = midString + std::to_string(OPP_main_loop_iter) + "_cell_acc.dat";

                opp_mpi_print_dat_to_txtfile(cell_e, f1.c_str());
                opp_mpi_print_dat_to_txtfile(cell_b, f2.c_str());
                opp_mpi_print_dat_to_txtfile(cell_j, f3.c_str());
                opp_mpi_print_dat_to_txtfile(cell_interp, f4.c_str());
                opp_mpi_print_dat_to_txtfile(cell_acc, f5.c_str());
#else           
                std::string f = std::string("F_") + midString + std::to_string(OPP_main_loop_iter);
                
                opp_print_dat_to_txtfile(cell_e, f.c_str(), "cell_e.dat");
                opp_print_dat_to_txtfile(cell_b, f.c_str(), "cell_b.dat");
                opp_print_dat_to_txtfile(cell_j, f.c_str(), "cell_j.dat");
                opp_print_dat_to_txtfile(cell_interp, f.c_str(), "cell_interp.dat");
                opp_print_dat_to_txtfile(cell_acc, f.c_str(), "cell_acc.dat");
#endif
            }

            // Leap frog method
            opp_loop_all__half_advance_b(
                cell_set,
                opp_get_arg(cell_ghost,                             OP_READ),
                opp_get_arg(cell_e, CellMap::xu_y_z, cell_cell_map, OP_WRITE),  // 0, adv_b_cell_map, 
                opp_get_arg(cell_e, CellMap::x_yu_z, cell_cell_map, OP_READ),   // 1, adv_b_cell_map, 
                opp_get_arg(cell_e, CellMap::x_y_zu, cell_cell_map, OP_READ),   // 2, adv_b_cell_map, 
                opp_get_arg(cell_e,                                 OP_READ),
                opp_get_arg(cell_b,                                 OP_INC)
            );

            serial_update_ghosts_B(cell_b);
            serial_update_ghosts(cell_j);
            serial_update_ghosts_B(cell_j);

            opp_loop_all__advance_e(
                cell_set,
                opp_get_arg(cell_iter_adv_e,                        OP_READ),
                opp_get_arg(cell_b, CellMap::xd_y_z, cell_cell_map, OP_READ), // 0, adv_e_cell_map
                opp_get_arg(cell_b, CellMap::x_yd_z, cell_cell_map, OP_READ), // 1, adv_e_cell_map
                opp_get_arg(cell_b, CellMap::x_y_zd, cell_cell_map, OP_READ), // 2, adv_e_cell_map
                opp_get_arg(cell_b,                                 OP_READ),
                opp_get_arg(cell_j,                                 OP_READ),
                opp_get_arg(cell_e,                                 OP_INC) 
            );

            opp_loop_all__half_advance_b(
                cell_set,
                opp_get_arg(cell_ghost,                             OP_READ),
                opp_get_arg(cell_e, CellMap::xu_y_z, cell_cell_map, OP_WRITE),  // 0, adv_b_cell_map, 
                opp_get_arg(cell_e, CellMap::x_yu_z, cell_cell_map, OP_READ),   // 1, adv_b_cell_map, 
                opp_get_arg(cell_e, CellMap::x_y_zu, cell_cell_map, OP_READ),   // 2, adv_b_cell_map, 
                opp_get_arg(cell_e,                                 OP_READ),
                opp_get_arg(cell_b,                                 OP_INC) 
            );

            serial_update_ghosts_B(cell_b);

            // if (opp_params->get<OPP_BOOL>("print_final"))
            // {
            //     log = get_global_level_log(node_charge_den, node_potential, particle_set->size, 
            //         n_parts_to_inject, (old_nparts - particle_set->size));
            // }
            
            if ((true && ((OPP_main_loop_iter+1) % 10 == 0 || OPP_main_loop_iter < 10)) ||
                OPP_main_loop_iter+1 == max_iter)
            {
                std::string f = std::string("F_") + std::to_string(OPP_main_loop_iter);
                
                opp_print_dat_to_txtfile(cell_e, f.c_str(), "cell_e.dat");
                opp_print_dat_to_txtfile(cell_b, f.c_str(), "cell_b.dat");
                opp_print_dat_to_txtfile(cell_j, f.c_str(), "cell_j.dat");

                opp_print_dat_to_txtfile(part_vel, f.c_str(), "part_vel.dat");
                opp_print_dat_to_txtfile(part_pos, f.c_str(), "part_pos.dat");
                opp_print_dat_to_txtfile(part_mesh_rel, f.c_str(), "part_mesh_rel.dat");
                opp_print_dat_to_txtfile(part_index, f.c_str(), "part_index.dat");

                std::string f1 = std::to_string(OPP_main_loop_iter) + "_cell_e.dat";
                std::string f2 = std::to_string(OPP_main_loop_iter) + "_cell_b.dat";
                std::string f3 = std::to_string(OPP_main_loop_iter) + "_cell_j.dat";
                std::string f4 = std::to_string(OPP_main_loop_iter) + "_part_vel.dat";
                std::string f5 = std::to_string(OPP_main_loop_iter) + "_part_pos.dat";
                std::string f6 = std::to_string(OPP_main_loop_iter) + "_part_mesh_rel.dat";

#ifdef USE_MPI
                opp_mpi_print_dat_to_txtfile(cell_e, f1.c_str());
                opp_mpi_print_dat_to_txtfile(cell_b, f2.c_str());
                opp_mpi_print_dat_to_txtfile(cell_j, f3.c_str());
                // opp_mpi_print_dat_to_txtfile(part_vel, f4.c_str());
                // opp_mpi_print_dat_to_txtfile(part_pos, f5.c_str());
                // opp_mpi_print_dat_to_txtfile(part_mesh_rel, f6.c_str());
#endif
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

