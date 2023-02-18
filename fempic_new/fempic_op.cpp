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

// AUTO GENERATED CODE

#include "fempic.h"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(-1);
    }

    opp::Params params(argv[1]);
    params.write(std::cout);

    FieldPointers mesh = LoadMesh(params);

    { // Start Scope for oppic
        oppic_init(argc, argv, &params);
 
        oppic_set nodes_set            = oppic_decl_set(mesh.n_nodes, "mesh_nodes");
        oppic_set cells_set            = oppic_decl_set(mesh.n_cells, "mesh_cells");
        oppic_set ifaces_set           = oppic_decl_set(mesh.n_ifaces, "inlet_faces_cells");
        oppic_set particles_set        = oppic_decl_particle_set("particles", cells_set); 

        oppic_map cell_to_nodes_map    = oppic_decl_map(cells_set,  nodes_set, NODES_PER_CELL,  mesh.cell_to_nodes,  "cell_to_nodes_map");
        oppic_map cell_to_cell_map     = oppic_decl_map(cells_set,  cells_set, NEIGHBOUR_CELLS, mesh.cell_to_cell,   "cell_to_cell_map"); 
        oppic_map iface_to_nodes_map   = oppic_decl_map(ifaces_set, nodes_set, 3,               mesh.iface_to_nodes, "iface_to_nodes_map");
        oppic_map iface_to_cell_map    = oppic_decl_map(ifaces_set, cells_set, 1,               mesh.iface_to_cell,  "iface_to_cell_map"); 

        oppic_dat cell_determinants    = oppic_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), "double", sizeof(double), (char*)mesh.cell_det,         "cell_determinants");  
        oppic_dat cell_volume          = oppic_decl_dat(cells_set, 1,                            "double", sizeof(double), (char*)mesh.cell_volume,      "cell_volume");        
        oppic_dat cell_electric_field  = oppic_decl_dat(cells_set, DIMENSIONS,                   "double", sizeof(double), (char*)mesh.cell_ef,          "cell_electric_field");
        oppic_dat cell_shape_deriv     = oppic_decl_dat(cells_set, (NODES_PER_CELL*DIMENSIONS),  "double", sizeof(double), (char*)mesh.cell_shape_deriv, "cell_shape_deriv"); 

        oppic_dat node_position        = oppic_decl_dat(nodes_set, DIMENSIONS, "double", sizeof(double), (char*)mesh.node_pos,     "node_position");      
        oppic_dat node_volume          = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)mesh.node_volume,  "node_volume");        
        oppic_dat node_bnd_potential   = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)mesh.node_bnd_pot, "node_bnd_potential"); 
        oppic_dat node_potential       = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)mesh.node_pot,     "node_potential");     
        oppic_dat node_charge_density  = oppic_decl_dat(nodes_set, 1,          "double", sizeof(double), (char*)mesh.node_ion_den, "node_charge_density");

        oppic_dat iface_v_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, "double", sizeof(double), (char*)mesh.iface_v_normal,      "iface_v_normal");        
        oppic_dat iface_u_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, "double", sizeof(double), (char*)mesh.iface_u_normal,      "iface_u_normal"); 
        oppic_dat iface_normal         = oppic_decl_dat(ifaces_set, DIMENSIONS, "double", sizeof(double), (char*)mesh.iface_normal,        "iface_normal");     
        oppic_dat iface_area           = oppic_decl_dat(ifaces_set, 1,          "double", sizeof(double), (char*)mesh.iface_area,          "iface_area");
        oppic_dat iface_inj_part_dist  = oppic_decl_dat(ifaces_set, 1,          "int",    sizeof(int),    (char*)mesh.iface_inj_part_dist, "iface_inj_part_dist");
        oppic_dat iface_node_pos       = oppic_decl_dat(ifaces_set, DIMENSIONS, "double", sizeof(double), (char*)mesh.iface_node_pos,      "iface_node_pos"); 

        oppic_dat part_position        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_position");
        oppic_dat part_velocity        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     "double", sizeof(double), nullptr, "part_velocity");    
        oppic_dat part_lc              = oppic_decl_particle_dat(particles_set, NODES_PER_CELL, "double", sizeof(double), nullptr, "part_lc");
        oppic_dat part_mesh_relation   = oppic_decl_particle_dat(particles_set, 1,              "int",    sizeof(int),    nullptr, "part_mesh_relation", true); // new cell index field

        mesh.DeleteValues();

        double remainder = 0.0;
        int max_iter = params.get<INT>("max_iter");
        // max_iter = 1;

        auto start = std::chrono::system_clock::now();

        for (int ts = 0; ts < max_iter; ts++)
        {
            int injected_count = 0;
            oppic_seq_loop_inject__Increase_particle_count(
                particles_set,                                                  // particles_set
                ifaces_set,                                                     // inlect_face_set
                oppic_arg_gbl(&(injected_count), 1, "int",   OP_RW),            // injected total global,
                oppic_arg_dat(iface_area,                    OP_READ),          // iface_area,
                oppic_arg_dat(iface_inj_part_dist,           OP_WRITE),         // iface_inj_part_dist,
                oppic_arg_gbl(&(remainder),     1, "double", OP_RW)             // remainder global,
            );

            int old_nparts = particles_set->size;
            oppic_par_loop_inject__InjectIons(
                particles_set,                                                              // particles_set
                oppic_arg_dat(part_position,                             OP_WRITE),         // part_position,
                oppic_arg_dat(part_velocity,                             OP_RW),            // part_velocity,
                oppic_arg_dat(part_mesh_relation,                        OP_WRITE),         // part_cell_connectivity,
                oppic_arg_dat(iface_to_cell_map,                         OP_READ, true),    // iface to cell map
                oppic_arg_dat(cell_electric_field, 0, iface_to_cell_map, OP_READ, true),    // cell_ef,
                oppic_arg_dat(iface_u_normal,                            OP_READ, true),    // iface_u,
                oppic_arg_dat(iface_v_normal,                            OP_READ, true),    // iface_v,
                oppic_arg_dat(iface_normal,                              OP_READ, true),    // iface_normal,
                oppic_arg_dat(iface_node_pos,                            OP_READ, true)     // iface_node_pos
            );

            oppic_reset_dat(node_charge_density, (char*)opp_zero_double16);
            oppic_par_loop_particle_all__MoveToCells(
                particles_set,                                                               // particles_set
                oppic_arg_dat(cell_electric_field,                       OP_READ, true),     // cell_ef,
                oppic_arg_dat(part_position,                             OP_RW),             // part_pos,
                oppic_arg_dat(part_velocity,                             OP_INC),            // part_vel,
                oppic_arg_dat(part_lc,                                   OP_RW),             // part_lc,
                oppic_arg_dat(part_mesh_relation,                        OP_RW),             // current_cell_index,
                oppic_arg_dat(cell_volume,                               OP_READ, true),     // current_cell_volume,
                oppic_arg_dat(cell_determinants,                         OP_READ, true),     // current_cell_det,
                oppic_arg_dat(cell_to_cell_map,                          OP_READ, true),     // cell_connectivity,
                oppic_arg_dat(node_charge_density, 0, cell_to_nodes_map, OP_INC,  true),     // node_charge_den0,
                oppic_arg_dat(node_charge_density, 1, cell_to_nodes_map, OP_INC,  true),     // node_charge_den1,
                oppic_arg_dat(node_charge_density, 2, cell_to_nodes_map, OP_INC,  true),     // node_charge_den2,
                oppic_arg_dat(node_charge_density, 3, cell_to_nodes_map, OP_INC,  true)      // node_charge_den3,
            );

            oppic_par_loop_all__ComputeNodeChargeDensity(
                nodes_set,                                       // nodes_set
                oppic_arg_dat(node_charge_density,  OP_RW),      // node_charge_density
                oppic_arg_dat(node_volume,          OP_READ)     // node_volume
            );

            // TODO: Add the field potential solver here

            oppic_reset_dat(cell_electric_field, (char*)opp_zero_double16);
            oppic_par_loop_all__ComputeElectricField(
                cells_set,                                                        // cells_set
                oppic_arg_dat(cell_electric_field,                  OP_INC),      // cell_electric_field,
                oppic_arg_dat(cell_shape_deriv,                     OP_READ),     // cell_shape_deriv,
                oppic_arg_dat(node_potential, 0, cell_to_nodes_map, OP_READ),     // node_potential0,
                oppic_arg_dat(node_potential, 1, cell_to_nodes_map, OP_READ),     // node_potential1,
                oppic_arg_dat(node_potential, 2, cell_to_nodes_map, OP_READ),     // node_potential2,
                oppic_arg_dat(node_potential, 3, cell_to_nodes_map, OP_READ)      // node_potential3,
            );
            
            {
                std::cout<<"ts: " << ts << "\t np: " << particles_set->size
                 <<" (" <<  injected_count << " added, "<< old_nparts - particles_set->size << " removed)" <<std::endl;
            }
        }

        std::chrono::duration<double> diff = std::chrono::system_clock::now() - start;
        std::cout << "\nFEMPIC - Time to iterate " << max_iter << " takes <chrono>: " << diff.count() << " s\n\n";

        oppic_exit();
    } // End Scope for oppic

    return 0;
}

//*************************************************************************************************

// {
//     std::string f = std::string("F_") + std::to_string(ts + 1);
//     oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
//     oppic_print_map_to_txtfile(cell_to_cell_map   , f.c_str(), "cell_to_cell_map.dat");
//     oppic_print_map_to_txtfile(iface_to_nodes_map , f.c_str(), "iface_to_nodes_map.dat");
//     oppic_print_map_to_txtfile(iface_to_cell_map  , f.c_str(), "iface_to_cell_map.dat");

//     oppic_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");
//     oppic_print_dat_to_txtfile(node_potential     , f.c_str(), "node_potential.dat");
//     oppic_print_dat_to_txtfile(cell_electric_field, f.c_str(), "cell_electric_field.dat");

//     oppic_print_dat_to_txtfile(iface_v_normal,      f.c_str(), "iface_v_normal.dat");
//     oppic_print_dat_to_txtfile(iface_u_normal,      f.c_str(), "iface_u_normal.dat");
//     oppic_print_dat_to_txtfile(iface_normal,        f.c_str(), "iface_normal.dat");
//     oppic_print_dat_to_txtfile(iface_area,          f.c_str(), "iface_area.dat");
//     oppic_print_dat_to_txtfile(iface_inj_part_dist, f.c_str(), "iface_inj_part_dist.dat");

//     oppic_print_dat_to_txtfile(part_position,      f.c_str(), "part_position.dat");
//     oppic_print_dat_to_txtfile(part_velocity,      f.c_str(), "part_velocity.dat");
//     oppic_print_dat_to_txtfile(part_lc,            f.c_str(), "part_lc.dat");
//     oppic_print_dat_to_txtfile(part_mesh_relation, f.c_str(), "part_mesh_relation.dat");
// }