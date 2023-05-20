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


#include "fempic.h"

void oppic_inject__Increase_particle_count(oppic_set,oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_inject__InjectIons(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_particle_all__MoveToCells(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);
void oppic_par_loop_all__ComputeNodeChargeDensity(oppic_set,oppic_arg,oppic_arg);
void oppic_par_loop_all__ComputeElectricField(oppic_set,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg,oppic_arg);

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
    //     exit(-1);
    // }

    opp::Params params("/ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_1.param"); //argv[1]);

    oppic_init(argc, argv, &params);

    params.write(std::cout);

{
    double plasma_den = params.get<OPP_REAL>("plasma_den");
    double dt = params.get<OPP_REAL>("dt");
    double ion_velocity = params.get<OPP_REAL>("ion_velocity");
    double spwt = params.get<OPP_REAL>("spwt");
    double mass = 2 * AMU;
    double charge = 1 * QE;
    int max_iter = params.get<OPP_INT>("max_iter");
    int ts = 0;
    int gn_nodes = 0, gn_cells = 0, gn_ifaces = 0; // global mesh counts

    std::shared_ptr<FieldPointers> g_mesh = std::make_shared<FieldPointers>();
    std::shared_ptr<FieldPointers> mesh = std::make_shared<FieldPointers>();
    if (OPP_my_rank == OPP_MPI_ROOT)
    {
        g_mesh = LoadMesh(params, argc, argv);
        gn_nodes  = g_mesh->n_nodes;
        gn_cells  = g_mesh->n_cells;
        gn_ifaces = g_mesh->n_ifaces;
    }

    MPI_Bcast(&gn_nodes, 1, MPI_INT, OPP_MPI_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&gn_cells, 1, MPI_INT, OPP_MPI_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&gn_ifaces, 1, MPI_INT, OPP_MPI_ROOT, MPI_COMM_WORLD);

    mesh->n_nodes  = opp_get_uniform_local_size(gn_nodes);
    mesh->n_cells  = opp_get_uniform_local_size(gn_cells);
    mesh->n_ifaces = opp_get_uniform_local_size(gn_ifaces);
    
    opp_printf("Main()", OPP_my_rank, "Global counts - Nodes[%d] Cells[%d] IFaces[%d]", gn_nodes, gn_cells, gn_ifaces);
    opp_printf("Main()", OPP_my_rank, "Before partitioning - Nodes[%d] Cells[%d] IFaces[%d]", mesh->n_nodes, mesh->n_cells, mesh->n_ifaces);

    mesh->CreateMeshArrays();

int* xxx = nullptr;

if (OPP_my_rank == OPP_MPI_ROOT)
{
    xxx = new int[gn_cells];
    for (int i = 0; i < gn_cells; i++) xxx[i] = i;
}
int* l_cells = new int[mesh->n_cells];
int* l_nodes = new int[mesh->n_nodes];
int* l_ifaces = new int[mesh->n_ifaces];

opp_uniform_scatter_array(xxx, l_cells, gn_cells, mesh->n_cells, 1); 
opp_uniform_scatter_array(xxx, l_nodes, gn_nodes, mesh->n_nodes, 1); 
opp_uniform_scatter_array(xxx, l_ifaces, gn_ifaces, mesh->n_ifaces, 1); 

    opp_uniform_scatter_array(g_mesh->cell_ef            , mesh->cell_ef            , gn_cells , mesh->n_cells , DIMENSIONS);
    opp_uniform_scatter_array(g_mesh->cell_to_nodes      , mesh->cell_to_nodes      , gn_cells , mesh->n_cells , NODES_PER_CELL); 
    opp_uniform_scatter_array(g_mesh->cell_to_cell       , mesh->cell_to_cell       , gn_cells , mesh->n_cells , NEIGHBOUR_CELLS); 
    opp_uniform_scatter_array(g_mesh->cell_det           , mesh->cell_det           , gn_cells , mesh->n_cells , DET_FIELDS * NEIGHBOUR_CELLS); 
    opp_uniform_scatter_array(g_mesh->cell_volume        , mesh->cell_volume        , gn_cells , mesh->n_cells , 1); 
    opp_uniform_scatter_array(g_mesh->cell_shape_deriv   , mesh->cell_shape_deriv   , gn_cells , mesh->n_cells , NODES_PER_CELL*DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->node_bnd_pot       , mesh->node_bnd_pot       , gn_nodes , mesh->n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_pot           , mesh->node_pot           , gn_nodes , mesh->n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_ion_den       , mesh->node_ion_den       , gn_nodes , mesh->n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_pos           , mesh->node_pos           , gn_nodes , mesh->n_nodes , DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->node_volume        , mesh->node_volume        , gn_nodes , mesh->n_nodes , 1);
    opp_uniform_scatter_array(g_mesh->node_type          , mesh->node_type          , gn_nodes , mesh->n_nodes , 1);
    opp_uniform_scatter_array(g_mesh->node_bnd_pot       , mesh->node_bnd_pot       , gn_nodes , mesh->n_nodes , 1);
    opp_uniform_scatter_array(g_mesh->iface_to_cell      , mesh->iface_to_cell      , gn_ifaces, mesh->n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_to_nodes     , mesh->iface_to_nodes     , gn_ifaces, mesh->n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_v_normal     , mesh->iface_v_normal     , gn_ifaces, mesh->n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_u_normal     , mesh->iface_u_normal     , gn_ifaces, mesh->n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_normal       , mesh->iface_normal       , gn_ifaces, mesh->n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_area         , mesh->iface_area         , gn_ifaces, mesh->n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_inj_part_dist, mesh->iface_inj_part_dist, gn_ifaces, mesh->n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_node_pos     , mesh->iface_node_pos     , gn_ifaces, mesh->n_ifaces, DIMENSIONS);       



    if (OPP_my_rank == OPP_MPI_ROOT)
    {
        g_mesh->DeleteValues();
    }

 
    oppic_set nodes_set            = oppic_decl_set(mesh->n_nodes, "mesh_nodes");
    oppic_set cells_set            = oppic_decl_set(mesh->n_cells, "mesh_cells");
    oppic_set ifaces_set           = oppic_decl_set(mesh->n_ifaces, "inlet_faces_cells");
    oppic_set particles_set        = oppic_decl_particle_set("particles", cells_set); 

    oppic_map cell_to_nodes_map    = oppic_decl_map(cells_set,  nodes_set, NODES_PER_CELL,  mesh->cell_to_nodes,  "cell_to_nodes_map");
    oppic_map cell_to_cell_map     = oppic_decl_map(cells_set,  cells_set, NEIGHBOUR_CELLS, mesh->cell_to_cell,   "cell_to_cell_map"); 
    oppic_map iface_to_cell_map    = oppic_decl_map(ifaces_set, cells_set, 1,               mesh->iface_to_cell,  "iface_to_cell_map"); 

// oppic_dat cell_to_cell_dat    = oppic_decl_dat(cells_set, NEIGHBOUR_CELLS, OPP_TYPE_INT, (char*)mesh->cell_to_cell, "cell_to_cell_dat");  
// oppic_dat cell_to_nodes_dat   = oppic_decl_dat(cells_set, NODES_PER_CELL, OPP_TYPE_INT, (char*)mesh->cell_to_nodes, "cell_to_nodes_dat");  
      
    oppic_dat cell_determinants    = oppic_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), OPP_TYPE_REAL, (char*)mesh->cell_det,         "cell_determinants");  
    oppic_dat cell_volume          = oppic_decl_dat(cells_set, 1,                            OPP_TYPE_REAL, (char*)mesh->cell_volume,      "cell_volume");        
    oppic_dat cell_electric_field  = oppic_decl_dat(cells_set, DIMENSIONS,                   OPP_TYPE_REAL, (char*)mesh->cell_ef,          "cell_electric_field");
    oppic_dat cell_shape_deriv     = oppic_decl_dat(cells_set, (NODES_PER_CELL*DIMENSIONS),  OPP_TYPE_REAL, (char*)mesh->cell_shape_deriv, "cell_shape_deriv"); 
oppic_dat cell_global_idx         = oppic_decl_dat(cells_set, 1,                            OPP_TYPE_INT, (char*)l_cells,      "cell_global_idx"); 
oppic_dat cell_part_count         = oppic_decl_dat(cells_set, 1,                            OPP_TYPE_INT, (char*)l_cells,      "cell_part_count"); // we just need to create a dat, reset before use

    oppic_dat node_volume          = oppic_decl_dat(nodes_set, 1,          OPP_TYPE_REAL, (char*)mesh->node_volume,  "node_volume");        
    oppic_dat node_potential       = oppic_decl_dat(nodes_set, 1,          OPP_TYPE_REAL, (char*)mesh->node_pot,     "node_potential");     
    oppic_dat node_charge_density  = oppic_decl_dat(nodes_set, 1,          OPP_TYPE_REAL, (char*)mesh->node_ion_den, "node_charge_density");
    oppic_dat node_pos             = oppic_decl_dat(nodes_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh->node_pos,     "node_pos");     
    oppic_dat node_type            = oppic_decl_dat(nodes_set, 1,          OPP_TYPE_INT,  (char*)mesh->node_type,    "node_type");
    oppic_dat node_bnd_pot         = oppic_decl_dat(nodes_set, 1,          OPP_TYPE_REAL, (char*)mesh->node_bnd_pot, "node_bnd_pot");
oppic_dat node_global_idx         = oppic_decl_dat(nodes_set, 1,                            OPP_TYPE_INT, (char*)l_nodes,      "node_global_idx"); 

    oppic_dat iface_v_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh->iface_v_normal,      "iface_v_normal");        
    oppic_dat iface_u_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh->iface_u_normal,      "iface_u_normal"); 
    oppic_dat iface_normal         = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh->iface_normal,        "iface_normal");     
    oppic_dat iface_area           = oppic_decl_dat(ifaces_set, 1,          OPP_TYPE_REAL, (char*)mesh->iface_area,          "iface_area");
    oppic_dat iface_inj_part_dist  = oppic_decl_dat(ifaces_set, 1,          OPP_TYPE_INT,  (char*)mesh->iface_inj_part_dist, "iface_inj_part_dist");
    oppic_dat iface_node_pos       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh->iface_node_pos,      "iface_node_pos"); 
oppic_dat iface_global_idx         = oppic_decl_dat(ifaces_set, 1,                            OPP_TYPE_INT, (char*)l_ifaces,      "iface_global_idx"); 

    oppic_dat part_position        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     OPP_TYPE_REAL, nullptr, "part_position");
    oppic_dat part_velocity        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     OPP_TYPE_REAL, nullptr, "part_velocity");    
    oppic_dat part_lc              = oppic_decl_particle_dat(particles_set, NODES_PER_CELL, OPP_TYPE_REAL, nullptr, "part_lc");
    oppic_dat part_mesh_relation   = oppic_decl_particle_dat(particles_set, 1,              OPP_TYPE_INT,  nullptr, "part_mesh_relation", true); // new cell index field

    oppic_set dummy_part_set       = oppic_decl_particle_set("dummy particles", cells_set); 
    oppic_dat dummy_part_random    = oppic_decl_particle_dat(dummy_part_set, 2, OPP_TYPE_REAL, nullptr, "dum_part_random");

    oppic_decl_const<OPP_REAL>(1, &spwt,         "CONST_spwt");
    oppic_decl_const<OPP_REAL>(1, &ion_velocity, "CONST_ion_velocity");
    oppic_decl_const<OPP_REAL>(1, &dt,           "CONST_dt");
    oppic_decl_const<OPP_REAL>(1, &plasma_den,   "CONST_plasma_den");
    oppic_decl_const<OPP_REAL>(1, &mass,         "CONST_mass");
    oppic_decl_const<OPP_REAL>(1, &charge,       "CONST_charge");

if (false)
{
    std::string f = std::string("Init_") + std::to_string(OPP_my_rank);
    oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
    // oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
    // oppic_print_map_to_txtfile(cell_to_cell_map  , f.c_str(), "cell_to_cell_map.dat");
    // oppic_print_dat_to_txtfile(cell_volume  , f.c_str(), "cell_volume.dat");
    // // oppic_print_dat_to_txtfile(cell_to_cell_dat  , f.c_str(), "cell_to_cell_dat.dat");
    // // oppic_print_dat_to_txtfile(cell_to_nodes_dat  , f.c_str(), "cell_to_nodes_dat.dat");
    // oppic_print_dat_to_txtfile(node_pos  , f.c_str(), "node_pos.dat");
    // oppic_print_dat_to_txtfile(node_type  , f.c_str(), "node_type.dat");
    // oppic_print_dat_to_txtfile(cell_shape_deriv  , f.c_str(), "cell_shape_deriv.dat");
    // oppic_print_dat_to_txtfile(iface_global_idx  , f.c_str(), "iface_global_idx.dat");
    // oppic_print_dat_to_txtfile(node_global_idx  , f.c_str(), "node_global_idx.dat");
    // oppic_print_dat_to_txtfile(cell_global_idx  , f.c_str(), "cell_global_idx.dat");
}

    mesh->DeleteValues();

    opp_partition(cells_set, cell_to_nodes_map);

    opp_printf("Main()", "After partitioning - Nodes[%d] Cells[%d] IFaces[%d]", mesh->n_nodes, mesh->n_cells, mesh->n_ifaces);

    int n_parts_to_inject = InitializeInjectDistributions(params, iface_inj_part_dist, iface_area, dummy_part_random);

    opp_printf("Main()", "After partitioning - Nodes[%d] Cells[%d] IFaces[%d] n_parts_to_inject[%d]", 
        nodes_set->size, cells_set->size, ifaces_set->size, n_parts_to_inject);

    std::shared_ptr<FESolver> solver = std::make_shared<FESolver>(&params, cell_to_nodes_map, node_type, node_pos, node_bnd_pot, argc, argv);
    solver->enrich_cell_shape_deriv(cell_shape_deriv);

if (false)
{
    std::string f = std::string("AfHalo_") + std::to_string(OPP_my_rank);
    // oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
    // oppic_print_map_to_txtfile(cell_to_cell_map  , f.c_str(), "cell_to_cell_map.dat");
    // oppic_print_dat_to_txtfile(cell_volume  , f.c_str(), "cell_volume.dat");
    // // oppic_print_dat_to_txtfile(cell_to_cell_dat  , f.c_str(), "cell_to_cell_dat.dat");
    // // oppic_print_dat_to_txtfile(cell_to_nodes_dat  , f.c_str(), "cell_to_nodes_dat.dat");
    // oppic_print_dat_to_txtfile(node_type  , f.c_str(), "node_type.dat");
    // oppic_print_dat_to_txtfile(node_bnd_pot  , f.c_str(), "node_bnd_pot.dat");
    oppic_print_dat_to_txtfile(dummy_part_random  , f.c_str(), "dummy_part_random.dat");
    oppic_print_dat_to_txtfile(iface_global_idx  , f.c_str(), "iface_global_idx.dat");
    oppic_print_dat_to_txtfile(node_global_idx  , f.c_str(), "node_global_idx.dat");
    oppic_print_dat_to_txtfile(cell_global_idx  , f.c_str(), "cell_global_idx.dat");
    // oppic_print_dat_to_txtfile(cell_part_count  , f.c_str(), "cell_part_count.dat");
}

    auto start = std::chrono::system_clock::now();
    auto start_iter1 = std::chrono::system_clock::now();

// opp_mpi_print_dat_to_txtfile(node_pos ,  "A_node_pos.dat");
// opp_mpi_print_dat_to_txtfile(node_volume ,  "A_node_volume.dat");

    for (ts = 0; ts < 100; ts++)
    {
        if (ts == 1) start_iter1 = std::chrono::system_clock::now();

        // if (OPP_my_rank == OPP_MPI_ROOT) 
            //opp_printf("Main()", "Starting main loop iteration %d XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", ts);

        opp_inc_part_count_with_distribution(particles_set, n_parts_to_inject, iface_inj_part_dist);

// std::string q = "Q_" + std::to_string(OPP_my_rank) + "-" + std::to_string(ts);
// oppic_print_dat_to_txtfile(part_mesh_relation, q.c_str(), "part_mesh_relation.dat");

        int old_nparts = particles_set->size;
        oppic_par_loop_inject__InjectIons(
            particles_set,                                                                               // particles_set
            oppic_arg_dat(part_position,                             OP_WRITE),                          // part_position,
            oppic_arg_dat(part_velocity,                             OP_WRITE),                          // part_velocity,
            oppic_arg_dat(part_mesh_relation,                        OP_RW),                             // part_cell_connectivity,
            oppic_arg_dat(iface_to_cell_map,                         OP_READ, OPP_Map_from_Mesh_Rel),    // iface to cell map
            oppic_arg_dat(cell_electric_field, 0, iface_to_cell_map, OP_READ, OPP_Map_from_Mesh_Rel),    // cell_ef,
            oppic_arg_dat(iface_u_normal,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_u,
            oppic_arg_dat(iface_v_normal,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_v,
            oppic_arg_dat(iface_normal,                              OP_READ, OPP_Map_from_Mesh_Rel),    // iface_normal,
            oppic_arg_dat(iface_node_pos,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_node_pos
            oppic_arg_dat(dummy_part_random,                         OP_READ, OPP_Map_from_Inj_part)     // dum_part_random
        );
// opp_printf("Main", "solve 2");
// {std::string f = std::string("K_r") + std::to_string(OPP_my_rank) + "_I" + std::to_string(ts);
// oppic_print_dat_to_txtfile(part_mesh_relation  , f.c_str(), "part_mesh_relation.dat");}

// std::string f = std::string("F_r") + std::to_string(OPP_my_rank) + std::string("_ts") + std::to_string(ts + 1);
// oppic_print_dat_to_txtfile(part_position  , f.c_str(), "part_position.dat");
// oppic_print_dat_to_txtfile(part_mesh_relation, f.c_str(), "part_mesh_relation.dat");

        oppic_reset_dat(node_charge_density, (char*)opp_zero_double16);
        oppic_par_loop_particle_all__MoveToCells(
            particles_set,                                                                                // particles_set
            oppic_arg_dat(cell_electric_field,                       OP_READ, OPP_Map_from_Mesh_Rel),     // cell_ef,
            oppic_arg_dat(part_position,                             OP_RW),                              // part_pos,
            oppic_arg_dat(part_velocity,                             OP_RW),                              // part_vel,
            oppic_arg_dat(part_lc,                                   OP_RW),                              // part_lc,
            oppic_arg_dat(part_mesh_relation,                        OP_RW),                              // current_cell_index,
            oppic_arg_dat(cell_volume,                               OP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_volume,
            oppic_arg_dat(cell_determinants,                         OP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_det,
            oppic_arg_dat(cell_to_cell_map,                          OP_READ, OPP_Map_from_Mesh_Rel),     // cell_connectivity,
            oppic_arg_dat(node_charge_density, 0, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den0,
            oppic_arg_dat(node_charge_density, 1, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den1,
            oppic_arg_dat(node_charge_density, 2, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den2,
            oppic_arg_dat(node_charge_density, 3, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel)      // node_charge_den3,
        );

// oppic_print_dat_to_txtfile(part_position  , f.c_str(), "AM_part_position.dat");
// oppic_print_dat_to_txtfile(part_mesh_relation, f.c_str(), "AM_part_mesh_relation.dat");
// oppic_print_dat_to_txtfile(part_velocity, f.c_str(), "AM_part_velocity.dat");
// oppic_print_dat_to_txtfile(node_charge_density  , f.c_str(), "AM_node_charge_density.dat");

// {std::string f = std::string("F_AM_mpi_node_charge_density_") + std::to_string(ts + 1);  
// opp_mpi_print_dat_to_txtfile(node_charge_density, f.c_str());}

    // {
    //     std::string f = std::to_string(ts) + std::string("_B_node_charge_density.dat");
    //     opp_mpi_print_dat_to_txtfile(node_charge_density, f.c_str());
    // }

// opp_printf("Main", "solve 3");
        oppic_par_loop_all__ComputeNodeChargeDensity(
            nodes_set,                                       // nodes_set
            oppic_arg_dat(node_charge_density,  OP_RW),      // node_charge_density
            oppic_arg_dat(node_volume,          OP_READ)     // node_volume
        );
    
    // {
    //     std::string f = std::to_string(ts) + std::string("_A_node_charge_density.dat");
    //     opp_mpi_print_dat_to_txtfile(node_charge_density, f.c_str());
    // }
// {std::string f = std::string("F_mpi_node_charge_density_") + std::to_string(ts + 1);  
// opp_mpi_print_dat_to_txtfile(node_charge_density, f.c_str());}
// {std::string f = std::string("RR") + std::to_string(OPP_my_rank) + "_I" + std::to_string(ts);
// oppic_print_dat_to_txtfile(node_charge_density  , f.c_str(), "node_charge_density.dat");}

// opp_printf("Main", "solve 4");
        solver->computePhi(  // TODO: Change this to kernel calls
            oppic_arg_dat(node_potential,      OP_WRITE),
            oppic_arg_dat(node_charge_density, OP_READ),
            oppic_arg_dat(node_bnd_pot,        OP_READ)
        );
// {std::string f = std::string("F_mpi_node_potential_") + std::to_string(ts + 1);  
// opp_mpi_print_dat_to_txtfile(node_potential, f.c_str());}
// opp_printf("Main", "solve 5");
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

// {std::string f = std::string("F_mpi_cell_electric_field_") + std::to_string(ts + 1);  
// opp_mpi_print_dat_to_txtfile(cell_electric_field, f.c_str());}

// {std::string f = std::string("F_mpi_node_potential_") + std::to_string(ts + 1);  
// opp_mpi_print_dat_to_txtfile(node_potential, f.c_str());}

// {std::string f = std::string("F_mpi_node_charge_density_") + std::to_string(ts + 1);        
// opp_mpi_print_dat_to_txtfile(node_charge_density, f.c_str());}


// opp_printf("Main", "solve 6");
// {std::string f = std::string("K_r") + std::to_string(OPP_my_rank) + "_I" + std::to_string(ts);
// oppic_print_dat_to_txtfile(cell_electric_field  , f.c_str(), "cell_electric_field.dat");
// oppic_print_dat_to_txtfile(node_potential  , f.c_str(), "node_potential.dat");
// oppic_print_dat_to_txtfile(cell_shape_deriv  , f.c_str(), "cell_shape_deriv.dat");}
        
        // if (false)
        {
            double max_den = 0.0, max_phi = 0.0;
            double global_max_den = 0.0, global_max_phi = 0.0;
            int global_part_size = 0, global_inj_size = 0, global_removed = 0;
            if (params.get<OPP_BOOL>("check_max_values"))
            {
                for (int n = 0; n< nodes_set->size; n++) // ideally, need to copy data from device to host, but at this point host has correct data
                {
                    // if (n+1 == nodes_set->size)
                    //     opp_printf("Main()", "%d XXXXXX max den: %2.25lE\t max |phi|: %2.25lE **** %d", n, max_den, ((double*)node_charge_density->data)[n], nodes_set->size);

                    if (abs(((double*)node_charge_density->data)[n]) > max_den) max_den = abs(((double*)node_charge_density->data)[n]);
                    if (abs(((double*)node_potential->data)[n]) > max_phi) max_phi = abs(((double*)node_potential->data)[n]);   
                }
            }

            // opp_printf("Main()", "OWN MAX max den: %2.25lE\t max |phi|: %2.25lE ****", max_den, max_phi);

            MPI_Reduce(&max_den, &global_max_den, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_ROOT, MPI_COMM_WORLD);
            MPI_Reduce(&max_phi, &global_max_phi, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_ROOT, MPI_COMM_WORLD);
            MPI_Reduce(&(particles_set->size), &global_part_size, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT, MPI_COMM_WORLD);
            MPI_Reduce(&n_parts_to_inject, &global_inj_size, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT, MPI_COMM_WORLD);
            int local_removed = (old_nparts - particles_set->size);
            MPI_Reduce(&local_removed, &global_removed, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT, MPI_COMM_WORLD);

            // opp_printf("Main()", "ts: %d\t np: %d (%d added, %d removed)\t max den: %2.5lE\t max |phi|: %2.5lE\t ****", 
            //     ts, particles_set->size, n_parts_to_inject, (old_nparts - particles_set->size), max_den, max_phi);

            if (OPP_my_rank == OPP_MPI_ROOT)
                opp_printf("Main()", "ts: %d\t np: %d (%d added, %d removed)\t max den: %2.25lE\t max |phi|: %2.5lE\t ****", 
                    ts, global_part_size, global_inj_size, global_removed, global_max_den, global_max_phi);
        }
    }

    std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
    std::chrono::duration<double> diff_1 = std::chrono::system_clock::now() - start_iter1;
    
    double d_diff = (double)diff.count();
    double d_diff_1 = (double)diff_1.count();
    double global_diff, global_diff_1;

    MPI_Reduce(&d_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&d_diff_1, &global_diff_1, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_ROOT, MPI_COMM_WORLD);

    // opp_printf("Main()", "FEMPIC - Time to iterate ts= %d <sec>: %lf and (ts-1) <sec>: %lf\n",
    //     max_iter, (double)diff.count(), (double)diff_1.count());
    if (OPP_my_rank == OPP_MPI_ROOT)
        opp_printf("Main()", "FEMPIC - Time to iterate ts= %d <sec>: %lf and (ts-1) <sec>: %lf\n",
            max_iter, global_diff, global_diff_1);

    // validate the final particle counts in each cell
    oppic_reset_dat(cell_part_count, (char*)opp_zero_int16);
    for (int p = 0; p < particles_set->size; p++)
    {
        int cell_index    = ((int *)part_mesh_relation->data)[p];
        ((int *)cell_part_count->data)[cell_index] += 1;
    }

    if (false)
    {
        opp_mpi_print_dat_to_txtfile(node_charge_density, "node_charge_density.dat");
        opp_mpi_print_dat_to_txtfile(node_potential, "node_potential.dat");
        opp_mpi_print_dat_to_txtfile(cell_electric_field, "cell_electric_field.dat");
        opp_mpi_print_dat_to_txtfile(cell_shape_deriv, "cell_shape_deriv.dat");
        opp_mpi_print_dat_to_txtfile(cell_determinants, "cell_determinants.dat");
        opp_mpi_print_dat_to_txtfile(cell_part_count, "cell_part_count.dat");
    }
}

    oppic_exit();


    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
// oppic_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");