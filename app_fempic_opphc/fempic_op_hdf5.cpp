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

#include "opp_hdf5.h"
#include "fempic_defs.h"
#include "fempic_misc.cpp"
#include "FESolver.h"

using namespace opp;

void opp_loop_all__init_boundary_pot(opp_set,opp_arg,opp_arg);
void opp_loop_inject__inject_ions(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
    opp_arg,opp_arg,opp_arg,opp_arg);
// void opp_particle_move__move_fused(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
//     opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_particle_move__move(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__compute_node_charge_density(opp_set,opp_arg,opp_arg);
void opp_loop_all__compute_electric_field(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__calc_pos_vel(opp_set,opp_arg,opp_arg,opp_arg);
void opp_loop_all__deposit_charge_on_nodes(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__get_max_values(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);

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
        std::string file      = opp_params->get<OPP_STRING>("hdf_filename");
        bool print_final_log  = opp_params->get<OPP_BOOL>("print_final");

        opp_set node_set         = opp_decl_set_hdf5(file.c_str(), "mesh_nodes");
        opp_set cell_set         = opp_decl_set_hdf5(file.c_str(), "mesh_cells");
        opp_set iface_set        = opp_decl_set_hdf5(file.c_str(), "inlet_faces_cells");
        opp_set particle_set     = opp_decl_particle_set_hdf5(file.c_str(), "particles", cell_set); 

        opp_map cell_v_nodes_map = opp_decl_map_hdf5(cell_set,  node_set, N_PER_C,  file.c_str(), "c_v_n_map");
        opp_map cell_v_cell_map  = opp_decl_map_hdf5(cell_set,  cell_set, NEIGHB_C, file.c_str(), "c_v_c_map"); 
        opp_map iface_v_cell_map = opp_decl_map_hdf5(iface_set, cell_set, ONE,      file.c_str(), "if_v_c_map"); 
        opp_map iface_v_node_map = opp_decl_map_hdf5(iface_set, node_set, N_PER_IF, file.c_str(), "if_v_n_map");

        opp_dat cell_det         = opp_decl_dat_hdf5(cell_set, ALL_DET,     DT_REAL, file.c_str(), "c_det");  
        opp_dat cell_volume      = opp_decl_dat_hdf5(cell_set, ONE,         DT_REAL, file.c_str(), "c_volume");        
        opp_dat cell_ef          = opp_decl_dat_hdf5(cell_set, DIM,         DT_REAL, file.c_str(), "c_ef");
        opp_dat cell_shape_deriv = opp_decl_dat_hdf5(cell_set, N_PER_C*DIM, DT_REAL, file.c_str(), "c_shape_deri"); 
        opp_dat global_cell_id   = opp_decl_dat_hdf5(cell_set, ONE,         DT_INT,  file.c_str(), "c_gbl_id"); 
        opp_dat cell_colors      = opp_decl_dat_hdf5(cell_set, ONE,         DT_INT,  file.c_str(), "c_colors");
        opp_dat cell_centroids   = opp_decl_dat_hdf5(cell_set, DIM,         DT_REAL, file.c_str(), "c_centroids");

        opp_dat node_volume      = opp_decl_dat_hdf5(node_set, ONE, DT_REAL, file.c_str(), "n_vol");        
        opp_dat node_potential   = opp_decl_dat_hdf5(node_set, ONE, DT_REAL, file.c_str(), "n_potential");     
        opp_dat node_charge_den  = opp_decl_dat_hdf5(node_set, ONE, DT_REAL, file.c_str(), "n_charge_den");
        opp_dat node_pos         = opp_decl_dat_hdf5(node_set, DIM, DT_REAL, file.c_str(), "n_pos");     
        opp_dat node_type        = opp_decl_dat_hdf5(node_set, ONE, DT_INT,  file.c_str(), "n_type");
        opp_dat node_bnd_pot     = opp_decl_dat_hdf5(node_set, ONE, DT_REAL, file.c_str(), "n_bnd_pot");
        // opp_dat node_id          = opp_decl_dat_hdf5(node_set, ONE, DT_INT,  file.c_str(), "n_id"); 
        // opp_dat node_colors      = opp_decl_dat_hdf5(node_set, ONE, DT_INT,  file.c_str(), "n_colors");

        opp_dat iface_v_norm  = opp_decl_dat_hdf5(iface_set, DIM,          DT_REAL, file.c_str(), "iface_v_norm");        
        opp_dat iface_u_norm  = opp_decl_dat_hdf5(iface_set, DIM,          DT_REAL, file.c_str(), "iface_u_norm"); 
        opp_dat iface_norm    = opp_decl_dat_hdf5(iface_set, DIM,          DT_REAL, file.c_str(), "iface_norm");     
        opp_dat iface_area    = opp_decl_dat_hdf5(iface_set, ONE,          DT_REAL, file.c_str(), "iface_area");
        opp_dat iface_dist    = opp_decl_dat_hdf5(iface_set, ONE,          DT_INT,  file.c_str(), "iface_dist");
        opp_dat iface_n_pos   = opp_decl_dat_hdf5(iface_set, N_PER_IF*DIM, DT_REAL, file.c_str(), "iface_n_pos"); 
        // opp_dat iface_id      = opp_decl_dat_hdf5(iface_set, ONE, DT_INT, file.c_str(), "iface_id"); 

        opp_dat part_position = opp_decl_dat_hdf5(particle_set, DIM,     DT_REAL, file.c_str(), "part_position");
        opp_dat part_velocity = opp_decl_dat_hdf5(particle_set, DIM,     DT_REAL, file.c_str(), "part_velocity");    
        opp_dat part_lc       = opp_decl_dat_hdf5(particle_set, N_PER_C, DT_REAL, file.c_str(), "part_lc");
        opp_dat part_mesh_rel = opp_decl_dat_hdf5(particle_set, ONE,     DT_INT,  file.c_str(), "part_mesh_rel", true);
        // opp_dat part_id       = opp_decl_dat_hdf5(particle_set, ONE,     DT_INT,  file.c_str(), "part_id");

        opp_set dummy_part_set   = opp_decl_particle_set_hdf5(file.c_str(), "dummy particles", cell_set); 
        opp_dat dummy_part_rand  = opp_decl_dat_hdf5(dummy_part_set, 2, DT_REAL, file.c_str(), "dummy_part_rand");

        opp_decl_const<OPP_REAL>(ONE, &spwt,           "CONST_spwt");
        opp_decl_const<OPP_REAL>(ONE, &ion_velocity,   "CONST_ion_velocity");
        opp_decl_const<OPP_REAL>(ONE, &dt,             "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &plasma_den,     "CONST_plasma_den");
        opp_decl_const<OPP_REAL>(ONE, &mass,           "CONST_mass");
        opp_decl_const<OPP_REAL>(ONE, &charge,         "CONST_charge");
        opp_decl_const<OPP_REAL>(ONE, &wall_potential, "CONST_wall_potential");

#ifdef USE_MPI
        genColoursForBlockPartition(cell_colors, cell_centroids, iface_n_pos, iface_v_node_map);

        // opp_partition(std::string("PARMETIS_KWAY"), cell_set, cell_v_nodes_map);
        // opp_partition(std::string("PARMETIS_GEOM"), iface_set, nullptr, iface_n_pos);
        // opp_partition(std::string("EXTERNAL"), node_set, nullptr, node_colors);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif
        
        opp_loop_all__init_boundary_pot(
            node_set, 
            opp_arg_dat(node_type, OPP_READ), 
            opp_arg_dat(node_bnd_pot, OPP_WRITE)
        );

        int n_parts_to_inject = InitializeInjectDistributions(iface_dist, iface_area, dummy_part_rand);

        initializeParticleMover(grid_spacing, DIM, node_pos, cell_volume, cell_det, global_cell_id);

        std::unique_ptr<FESolver> field_solver = std::make_unique<FESolver>(cell_v_nodes_map, 
            node_type, node_pos, node_bnd_pot, argc, argv);
            
        field_solver->enrich_cell_shape_deriv(cell_shape_deriv);

        opp_inc_part_count_with_distribution(particle_set, n_parts_to_inject, iface_dist, false);

        opp_profiler->end("Setup");

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        int64_t total_part_iter = 0;
        opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            if (OPP_DBG && OPP_rank == OPP_ROOT) 
                opp_printf("Main", "Starting main loop iteration %d *************", OPP_main_loop_iter);

            if (OPP_main_loop_iter != 0)
                opp_inc_part_count_with_distribution(particle_set, n_parts_to_inject, iface_dist, false);

            int old_nparts = particle_set->size;
            opp_loop_inject__inject_ions(
                particle_set,                                                                           
                opp_arg_dat(part_position,                OPP_WRITE),                      
                opp_arg_dat(part_velocity,                OPP_WRITE),                      
                opp_arg_dat(part_mesh_rel,                OPP_RW),
                opp_arg_dat(part_mesh_rel,                OPP_RW),
                opp_arg_dat(iface_v_cell_map,             OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(cell_ef, 0, iface_v_cell_map, OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(iface_u_norm,                 OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(iface_v_norm,                 OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(iface_norm,                   OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(iface_n_pos,                  OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(dummy_part_rand,              OPP_READ, OPP_Map_from_Inj_part) 
            );

            opp_reset_dat(
                node_charge_den, 
                (char*)opp_zero_double16);

            opp_loop_all__calc_pos_vel(
                particle_set,                                                                           
                opp_arg_dat(cell_ef,       OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(part_position, OPP_WRITE),                         
                opp_arg_dat(part_velocity, OPP_WRITE)
            );

            opp_particle_move__move(
                particle_set,
                opp_arg_dat(part_position,   OPP_READ),
                opp_arg_dat(part_mesh_rel,   OPP_RW),
                opp_arg_dat(part_lc,         OPP_WRITE),
                opp_arg_dat(cell_volume,     OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(cell_det,        OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(cell_v_cell_map, OPP_READ, OPP_Map_from_Mesh_Rel)
            );

            opp_loop_all__deposit_charge_on_nodes(
                particle_set, 
                opp_arg_dat(part_lc,                              OPP_READ),
                opp_arg_dat(node_charge_den, 0, cell_v_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),
                opp_arg_dat(node_charge_den, 1, cell_v_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),
                opp_arg_dat(node_charge_den, 2, cell_v_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),
                opp_arg_dat(node_charge_den, 3, cell_v_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel)
            );

            opp_loop_all__compute_node_charge_density(
                node_set,                            
                opp_arg_dat(node_charge_den,  OPP_RW), 
                opp_arg_dat(node_volume,      OPP_READ)
            );

            field_solver->computePhi(  // TODO: Change this to kernel calls
                opp_arg_dat(node_potential,  OPP_WRITE),
                opp_arg_dat(node_charge_den, OPP_READ),
                opp_arg_dat(node_bnd_pot,    OPP_READ)
            );

            opp_reset_dat(
                cell_ef, 
                (char*)opp_zero_double16); 

            opp_loop_all__compute_electric_field(
                cell_set,                                                   
                opp_arg_dat(cell_ef,                             OPP_INC), 
                opp_arg_dat(cell_shape_deriv,                    OPP_READ),
                opp_arg_dat(node_potential, 0, cell_v_nodes_map, OPP_READ),
                opp_arg_dat(node_potential, 1, cell_v_nodes_map, OPP_READ),
                opp_arg_dat(node_potential, 2, cell_v_nodes_map, OPP_READ),
                opp_arg_dat(node_potential, 3, cell_v_nodes_map, OPP_READ) 
            );

            if (print_final_log || OPP_main_loop_iter + 1 == max_iter)
            {
                OPP_REAL max_n_chg_den = 0.0, max_n_pot = 0.0;

                opp_loop_all__get_max_values(
                    node_set,
                    opp_arg_dat(node_charge_den, OPP_READ),
                    opp_arg_gbl(&max_n_chg_den, 1, "double", OPP_MAX),
                    opp_arg_dat(node_potential, OPP_READ),
                    opp_arg_gbl(&max_n_pot, 1, "double", OPP_MAX)
                );

                log = get_global_level_log(max_n_chg_den, max_n_pot, particle_set->size, 
                    n_parts_to_inject, (old_nparts - particle_set->size));
            }
            // opp_printf("XXX", "particle_set->size=%d", particle_set->size); 
            
            total_part_iter += particle_set->size;  

            if ((print_final_log || OPP_main_loop_iter + 1 == max_iter) && OPP_rank == OPP_ROOT) 
                opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
        opp_profiler->end("MainLoop");

        // if (OPP_DBG)
        //     print_per_cell_particle_counts(cell_colors, part_mesh_rel); // cell_colors will reset

        opp_printf("Main","total particles= %" PRId64, total_part_iter); 

        int64_t gbl_total_part_iter = 0;
#ifdef USE_MPI
        MPI_Reduce(&total_part_iter, &gbl_total_part_iter, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
#else
        gbl_total_part_iter = total_part_iter;
#endif
        if (OPP_rank == OPP_ROOT) 
            opp_printf("Main", "Main loop completed after %d iterations with %" PRId64 " particle iterations ****", 
                max_iter, gbl_total_part_iter);  
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(cell_v_nodes_map  , f.c_str(), "cell_v_nodes_map.dat");
// opp_print_dat_to_txtfile(node_charge_den, f.c_str(), "node_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(cell_shape_deriv, "cell_shape_deriv.dat");