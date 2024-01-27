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

#include "fempic_loader.h"
#include "fempic_misc.cpp"
#include "FESolver.h"

using namespace opp;

void opp_loop_all__InitBndPotential(opp_set,opp_arg,opp_arg);
void opp_loop_inject__InjectIons(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
    opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all_part_move__MoveToCells(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
    opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__ComputeNodeChargeDensity(opp_set,opp_arg,opp_arg);
void opp_loop_all__ComputeElectricField(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__CalculateNewPartPosVel(opp_set,opp_arg,opp_arg,opp_arg);
void opp_loop_all__DepositChargeOnNodes(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_particle_mover__Move(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);

void initializeParticleMover(const double gridSpacing, int dim, const opp_dat node_pos_dat, 
    const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat global_cell_id_dat);

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << 
        "/ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_1.param" << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    {
        opp_profiler->start("Setup");

        double plasma_den     = opp_params->get<OPP_REAL>("plasma_den");
        double dt             = opp_params->get<OPP_REAL>("dt");
        double ion_velocity   = opp_params->get<OPP_REAL>("ion_velocity");
        double spwt           = opp_params->get<OPP_REAL>("spwt");
        double wall_potential = opp_params->get<OPP_REAL>("wall_potential");
        double grid_spacing   = opp_params->get<OPP_REAL>("grid_spacing");
        double mass           = 2 * AMU;
        double charge         = 1 * QE;
        int max_iter          = opp_params->get<OPP_INT>("max_iter");   
        std::string log       = "";

        std::shared_ptr<FieldPointers> g_m, m; // g_m - global mesh, m - local mesh
        g_m = std::make_shared<FieldPointers>();
        
        if (OPP_rank == OPP_ROOT)
        {
            // Load using the original FemPIC loaders and distribute
            g_m = LoadMesh(opp_params->get<OPP_STRING>("global_mesh"), opp_params->get<OPP_STRING>("inlet_mesh"),
                                opp_params->get<OPP_STRING>("wall_mesh"), opp_params->get<OPP_STRING>("cluster"));
            
            opp_printf("Main", "Global counts - Nodes[%d] Cells[%d] IFaces[%d]", 
                g_m->n_nodes, g_m->n_cells, g_m->n_ifaces);
        }

        DistributeMeshOverRanks(g_m, m);

        opp_set node_set         = opp_decl_mesh_set(m->n_nodes, "mesh_nodes");
        opp_set cell_set         = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_set iface_set        = opp_decl_mesh_set(m->n_ifaces, "inlet_faces_cells");
        opp_set particle_set     = opp_decl_part_set("particles", cell_set); 

        opp_map cell_v_nodes_map = opp_decl_mesh_map(cell_set,  node_set, N_PER_C,  m->c_to_n, "c_v_n_map");
        opp_map cell_v_cell_map  = opp_decl_mesh_map(cell_set,  cell_set, NEIGHB_C, m->c_to_c,  "c_v_c_map"); 
        opp_map iface_v_cell_map = opp_decl_mesh_map(iface_set, cell_set, ONE,      m->if_to_c, "if_v_c_map"); 
        opp_map iface_v_node_map = opp_decl_mesh_map(iface_set, node_set, N_PER_IF, m->if_to_n, "if_v_n_map");

        opp_dat cell_det         = opp_decl_mesh_dat(cell_set, ALL_DET,     DT_REAL, m->c_det,      "c_det");  
        opp_dat cell_volume      = opp_decl_mesh_dat(cell_set, ONE,         DT_REAL, m->c_vol,      "c_volume");        
        opp_dat cell_ef          = opp_decl_mesh_dat(cell_set, DIM,         DT_REAL, m->c_ef,       "c_ef");
        opp_dat cell_shape_deriv = opp_decl_mesh_dat(cell_set, N_PER_C*DIM, DT_REAL, m->c_sd,       "c_shape_deri"); 
        opp_dat global_cell_id   = opp_decl_mesh_dat(cell_set, ONE,         DT_INT,  m->c_id,       "c_gbl_id"); 
        opp_dat cell_colors      = opp_decl_mesh_dat(cell_set, ONE,         DT_INT,  m->c_col,      "c_colors");
        opp_dat cell_centroids   = opp_decl_mesh_dat(cell_set, DIM,         DT_REAL, m->c_centroid, "c_centroids");

        opp_dat node_volume      = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_vol,     "n_vol");        
        opp_dat node_potential   = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_pot,     "n_potential");     
        opp_dat node_charge_den  = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_ion_den, "n_charge_den");
        opp_dat node_pos         = opp_decl_mesh_dat(node_set, DIM, DT_REAL, m->n_pos,     "n_pos");     
        opp_dat node_type        = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_type,    "n_type");
        opp_dat node_bnd_pot     = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_bnd_pot, "n_bnd_pot");
        // opp_dat node_id          = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_id,      "n_id"); 
        // opp_dat node_colors      = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_color,   "n_colors");

        opp_dat iface_v_norm  = opp_decl_mesh_dat(iface_set, DIM,          DT_REAL, m->if_v_norm, "iface_v_norm");        
        opp_dat iface_u_norm  = opp_decl_mesh_dat(iface_set, DIM,          DT_REAL, m->if_u_norm, "iface_u_norm"); 
        opp_dat iface_norm    = opp_decl_mesh_dat(iface_set, DIM,          DT_REAL, m->if_norm,   "iface_norm");     
        opp_dat iface_area    = opp_decl_mesh_dat(iface_set, ONE,          DT_REAL, m->if_area,   "iface_area");
        opp_dat iface_dist    = opp_decl_mesh_dat(iface_set, ONE,          DT_INT,  m->if_dist,   "iface_dist");
        opp_dat iface_n_pos   = opp_decl_mesh_dat(iface_set, N_PER_IF*DIM, DT_REAL, m->if_n_pos,  "iface_n_pos"); 
        // opp_dat iface_id      = opp_decl_mesh_dat(iface_set, ONE, DT_INT,  m->if_id,     "iface_id"); 

        opp_dat part_position = opp_decl_part_dat(particle_set, DIM,     DT_REAL, nullptr, "part_position");
        opp_dat part_velocity = opp_decl_part_dat(particle_set, DIM,     DT_REAL, nullptr, "part_velocity");    
        opp_dat part_lc       = opp_decl_part_dat(particle_set, N_PER_C, DT_REAL, nullptr, "part_lc");
        opp_dat part_mesh_rel = opp_decl_part_dat(particle_set, ONE,     DT_INT,  nullptr, "part_mesh_rel", true);
        // opp_dat part_id       = opp_decl_part_dat(particle_set, ONE,     DT_INT,  nullptr, "part_id");

        opp_set dummy_part_set   = opp_decl_part_set("dummy particles", cell_set); 
        opp_dat dummy_part_rand  = opp_decl_part_dat(dummy_part_set, 2, DT_REAL, nullptr, "dummy_part_rand");

        opp_decl_const<OPP_REAL>(ONE, &spwt,           "CONST_spwt");
        opp_decl_const<OPP_REAL>(ONE, &ion_velocity,   "CONST_ion_velocity");
        opp_decl_const<OPP_REAL>(ONE, &dt,             "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &plasma_den,     "CONST_plasma_den");
        opp_decl_const<OPP_REAL>(ONE, &mass,           "CONST_mass");
        opp_decl_const<OPP_REAL>(ONE, &charge,         "CONST_charge");
        opp_decl_const<OPP_REAL>(ONE, &wall_potential, "CONST_wall_potential");

        m->DeleteValues();

#ifdef USE_MPI
        genColoursForBlockPartition(cell_colors, cell_centroids, iface_n_pos, iface_v_node_map);

        // opp_partition(std::string("PARMETIS_KWAY"), cell_set, cell_v_nodes_map);
        // opp_partition(std::string("PARMETIS_GEOM"), iface_set, nullptr, iface_n_pos);
        // opp_partition(std::string("EXTERNAL"), node_set, nullptr, node_colors);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif
        
        opp_loop_all__InitBndPotential(
            node_set, 
            opp_get_arg(node_type, OP_READ), 
            opp_get_arg(node_bnd_pot, OP_WRITE)
        );

        int n_parts_to_inject = InitializeInjectDistributions(iface_dist, iface_area, dummy_part_rand);

        initializeParticleMover(grid_spacing, DIM, node_pos, cell_volume, cell_det, global_cell_id);

        std::unique_ptr<FESolver> field_solver = std::make_unique<FESolver>(cell_v_nodes_map, 
            node_type, node_pos, node_bnd_pot, argc, argv);
            
        field_solver->enrich_cell_shape_deriv(cell_shape_deriv);

        opp_profiler->end("Setup");

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif


        opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            if (OP_DEBUG && OPP_rank == OPP_ROOT) 
                opp_printf("Main", "Starting main loop iteration %d *************", OPP_main_loop_iter);

            opp_inc_part_count_with_distribution(particle_set, n_parts_to_inject, iface_dist, false);

            int old_nparts = particle_set->size;
            opp_loop_inject__InjectIons(
                particle_set,                                                                           
                opp_get_arg(part_position,                OP_WRITE),                      
                opp_get_arg(part_velocity,                OP_WRITE),                      
                opp_get_arg(part_mesh_rel,                OP_RW),
                opp_get_arg(part_mesh_rel,                OP_RW),
                opp_get_arg(iface_v_cell_map,             OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_ef, 0, iface_v_cell_map, OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(iface_u_norm,                 OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(iface_v_norm,                 OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(iface_norm,                   OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(iface_n_pos,                  OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(dummy_part_rand,              OP_READ, OPP_Map_from_Inj_part) 
            );

            opp_reset_dat(
                node_charge_den, 
                (char*)opp_zero_double16);

            opp_loop_all__CalculateNewPartPosVel(
                particle_set,                                                                           
                opp_get_arg(cell_ef,       OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(part_position, OP_WRITE),                         
                opp_get_arg(part_velocity, OP_WRITE)
            );

            opp_particle_mover__Move(
                particle_set,
                opp_get_arg(part_position,   OP_READ),
                opp_get_arg(part_mesh_rel,   OP_RW),
                opp_get_arg(part_lc,         OP_WRITE),
                opp_get_arg(cell_volume,     OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_det,        OP_READ, OPP_Map_from_Mesh_Rel),
                opp_get_arg(cell_v_cell_map, OP_READ, OPP_Map_from_Mesh_Rel)
            );

            opp_loop_all__DepositChargeOnNodes(
                particle_set, 
                opp_get_arg(part_lc,                              OP_READ),
                opp_get_arg(node_charge_den, 0, cell_v_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),
                opp_get_arg(node_charge_den, 1, cell_v_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),
                opp_get_arg(node_charge_den, 2, cell_v_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),
                opp_get_arg(node_charge_den, 3, cell_v_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel)
            );

            opp_loop_all__ComputeNodeChargeDensity(
                node_set,                            
                opp_get_arg(node_charge_den,  OP_RW), 
                opp_get_arg(node_volume,      OP_READ)
            );

            field_solver->computePhi(  // TODO: Change this to kernel calls
                opp_get_arg(node_potential,  OP_WRITE),
                opp_get_arg(node_charge_den, OP_READ),
                opp_get_arg(node_bnd_pot,    OP_READ)
            );

            opp_reset_dat(
                cell_ef, 
                (char*)opp_zero_double16); 

            opp_loop_all__ComputeElectricField(
                cell_set,                                                   
                opp_get_arg(cell_ef,                             OP_INC), 
                opp_get_arg(cell_shape_deriv,                    OP_READ),
                opp_get_arg(node_potential, 0, cell_v_nodes_map, OP_READ),
                opp_get_arg(node_potential, 1, cell_v_nodes_map, OP_READ),
                opp_get_arg(node_potential, 2, cell_v_nodes_map, OP_READ),
                opp_get_arg(node_potential, 3, cell_v_nodes_map, OP_READ) 
            );

            if (opp_params->get<OPP_BOOL>("print_final"))
            {
                log = get_global_level_log(node_charge_den, node_potential, particle_set->size, 
                    n_parts_to_inject, (old_nparts - particle_set->size));
            }
            
            if (OPP_rank == OPP_ROOT)
                opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
        opp_profiler->end("MainLoop");

        if (OP_DEBUG)
            print_per_cell_particle_counts(cell_colors, part_mesh_rel); // cell_colors will reset
        
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Main loop completed after %d iterations ****", max_iter);
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(cell_v_nodes_map  , f.c_str(), "cell_v_nodes_map.dat");
// opp_print_dat_to_txtfile(node_charge_den, f.c_str(), "node_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(cell_shape_deriv, "cell_shape_deriv.dat");

// std::string f = std::string("F_") + std::to_string(OPP_rank);
// if (OPP_rank < 10) opp_print_dat_to_txtfile(cell_colors, f.c_str(), "cell_colors_BEFORE.dat");
// opp_reset_dat(cell_colors, (char*)opp_zero_int16);
// genColoursForBlockPartition(cell_colors, cell_centroids, iface_n_pos, iface_v_node_map);
// if (OPP_rank < 10) opp_print_dat_to_txtfile(cell_colors, f.c_str(), "cell_colors_AFTER.dat");