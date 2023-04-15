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
#include <opp_mpi.h>

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

    opp::Params params("/ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_1.param");
    
    oppic_init(argc, argv, &params);
    
    opp_printf("Main()", "RUNNING fempic TEST");  

    params.write(std::cout);

    std::shared_ptr<FieldPointers> g_mesh, g1_mesh;

    g1_mesh = LoadMesh(params, argc, argv);

opp_printf("XXXXXXXXXX", OPP_my_rank, " 00");  
    
    g1_mesh->DeleteValues();

opp_printf("XXXXXXXXXX", OPP_my_rank, " 1");

    if (OPP_my_rank == OPP_MPI_ROOT)
    {
        g_mesh = LoadMesh(params, argc, argv);
    }
    else
    {
        g_mesh = std::make_shared<FieldPointers>();
    }

opp_printf("XXXXXXXXXX", OPP_my_rank, " 2");
    
    FieldPointers mesh;

    mesh.n_nodes  = opp_get_uniform_local_size(g1_mesh->n_nodes);
    mesh.n_cells  = opp_get_uniform_local_size(g1_mesh->n_cells);
    mesh.n_ifaces = opp_get_uniform_local_size(g1_mesh->n_ifaces);

opp_printf("opp_get_uniform_local_size", OPP_my_rank, " %d %d %d",mesh.n_nodes,mesh.n_cells,mesh.n_ifaces);
opp_printf("XXXXXXXXXX", OPP_my_rank, " 3");

    mesh.CreateMeshArrays();

opp_printf("XXXXXXXXXX", OPP_my_rank, " 4");

    opp_uniform_scatter_array(g_mesh->cell_ef            , mesh.cell_ef            , g1_mesh->n_cells , mesh.n_cells , DIMENSIONS);
    opp_uniform_scatter_array(g_mesh->cell_to_nodes      , mesh.cell_to_nodes      , g1_mesh->n_cells , mesh.n_cells , NODES_PER_CELL); 
    opp_uniform_scatter_array(g_mesh->cell_to_cell       , mesh.cell_to_cell       , g1_mesh->n_cells , mesh.n_cells , NEIGHBOUR_CELLS); 
    opp_uniform_scatter_array(g_mesh->cell_det           , mesh.cell_det           , g1_mesh->n_cells , mesh.n_cells , DET_FIELDS * NEIGHBOUR_CELLS); 
    opp_uniform_scatter_array(g_mesh->cell_volume        , mesh.cell_volume        , g1_mesh->n_cells , mesh.n_cells , 1); 
    opp_uniform_scatter_array(g_mesh->cell_shape_deriv   , mesh.cell_shape_deriv   , g1_mesh->n_cells , mesh.n_cells , NODES_PER_CELL*DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->node_bnd_pot       , mesh.node_bnd_pot       , g1_mesh->n_nodes , mesh.n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_pot           , mesh.node_pot           , g1_mesh->n_nodes , mesh.n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_ion_den       , mesh.node_ion_den       , g1_mesh->n_nodes , mesh.n_nodes , 1); 
    opp_uniform_scatter_array(g_mesh->node_pos           , mesh.node_pos           , g1_mesh->n_nodes , mesh.n_nodes , DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->node_volume        , mesh.node_volume        , g1_mesh->n_nodes , mesh.n_nodes , 1);
    opp_uniform_scatter_array(g_mesh->iface_to_cell      , mesh.iface_to_cell      , g1_mesh->n_ifaces, mesh.n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_to_nodes     , mesh.iface_to_nodes     , g1_mesh->n_ifaces, mesh.n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_v_normal     , mesh.iface_v_normal     , g1_mesh->n_ifaces, mesh.n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_u_normal     , mesh.iface_u_normal     , g1_mesh->n_ifaces, mesh.n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_normal       , mesh.iface_normal       , g1_mesh->n_ifaces, mesh.n_ifaces, DIMENSIONS); 
    opp_uniform_scatter_array(g_mesh->iface_area         , mesh.iface_area         , g1_mesh->n_ifaces, mesh.n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_inj_part_dist, mesh.iface_inj_part_dist, g1_mesh->n_ifaces, mesh.n_ifaces, 1); 
    opp_uniform_scatter_array(g_mesh->iface_node_pos     , mesh.iface_node_pos     , g1_mesh->n_ifaces, mesh.n_ifaces, DIMENSIONS);        

opp_printf("XXXXXXXXXX", OPP_my_rank, " 5");

    if (OPP_my_rank == OPP_MPI_ROOT)
    {
        g_mesh->DeleteValues();
    }

opp_printf("XXXXXXXXXX", OPP_my_rank, " 6");

    double plasma_den = params.get<OPP_REAL>("plasma_den");
    double dt = params.get<OPP_REAL>("dt");
    double ion_velocity = params.get<OPP_REAL>("ion_velocity");
    double spwt = params.get<OPP_REAL>("spwt");
    double mass = 2 * AMU;
    double charge = 1 * QE;
    int max_iter = params.get<OPP_INT>("max_iter");
    double remainder = 0.0;
    int ts = 0;

opp_printf("main()", OPP_my_rank, " starting to decl opp structs: nodes:%d cells:%d, ifaces:%d", mesh.n_nodes, mesh.n_cells, mesh.n_ifaces);
//     { // Start Scope for oppic

        int pcount = 20;
        double* arr_part_position    = new double[pcount * DIMENSIONS];
        double* arr_part_velocity    = new double[pcount * DIMENSIONS];
        int* arr_part_lc             = new int[pcount * NODES_PER_CELL];
        int* arr_part_mesh_relation  = new int[pcount];

        for (int z = 0; z < pcount; z++)
        {
            arr_part_position[z * DIMENSIONS + 0] = (OPP_my_rank+1) * 1000 + z * 10 + 0;
            arr_part_position[z * DIMENSIONS + 1] = (OPP_my_rank+1) * 1000 + z * 10 + 1;
            arr_part_position[z * DIMENSIONS + 2] = (OPP_my_rank+1) * 1000 + z * 10 + 2;

            arr_part_velocity[z * DIMENSIONS + 0] = (OPP_my_rank+1) * 10000 + z * 10 + 0;
            arr_part_velocity[z * DIMENSIONS + 1] = (OPP_my_rank+1) * 10000 + z * 10 + 1;
            arr_part_velocity[z * DIMENSIONS + 2] = (OPP_my_rank+1) * 10000 + z * 10 + 2;

            arr_part_lc[z * NODES_PER_CELL + 0] = (OPP_my_rank+1) * 100000 + z * 10 + 0;
            arr_part_lc[z * NODES_PER_CELL + 1] = (OPP_my_rank+1) * 100000 + z * 10 + 1;
            arr_part_lc[z * NODES_PER_CELL + 2] = (OPP_my_rank+1) * 100000 + z * 10 + 2;
            arr_part_lc[z * NODES_PER_CELL + 3] = (OPP_my_rank+1) * 100000 + z * 10 + 3;

            arr_part_mesh_relation[z] = ((z*873 +1) % mesh.n_cells);
        }
 
        oppic_set nodes_set            = oppic_decl_set(mesh.n_nodes, "mesh_nodes");
        oppic_set cells_set            = oppic_decl_set(mesh.n_cells, "mesh_cells");
        oppic_set ifaces_set           = oppic_decl_set(mesh.n_ifaces, "mesh_inlet_faces");
        
// opp_printf("XXXXXXXXXX", OPP_my_rank, " 62");
        oppic_map cell_to_nodes_map    = oppic_decl_map(cells_set,  nodes_set, NODES_PER_CELL,  mesh.cell_to_nodes,  "cell_to_nodes_map");
        oppic_map cell_to_cell_map     = oppic_decl_map(cells_set,  cells_set, NEIGHBOUR_CELLS, mesh.cell_to_cell,   "cell_to_cell_map"); 
        oppic_map iface_to_cell_map    = oppic_decl_map(ifaces_set, cells_set, 1,               mesh.iface_to_cell,  "iface_to_cell_map"); 
// opp_printf("XXXXXXXXXX", OPP_my_rank, " 63");
        oppic_dat cell_determinants    = oppic_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), OPP_TYPE_REAL, (char*)mesh.cell_det,         "cell_determinants");  
        oppic_dat cell_volume          = oppic_decl_dat(cells_set, 1,                            OPP_TYPE_REAL, (char*)mesh.cell_volume,      "cell_volume");        
        oppic_dat cell_electric_field  = oppic_decl_dat(cells_set, DIMENSIONS,                   OPP_TYPE_REAL, (char*)mesh.cell_ef,          "cell_electric_field");
        oppic_dat cell_shape_deriv     = oppic_decl_dat(cells_set, (NODES_PER_CELL*DIMENSIONS),  OPP_TYPE_REAL, (char*)mesh.cell_shape_deriv, "cell_shape_deriv"); 
// opp_printf("XXXXXXXXXX", OPP_my_rank, " 64");   
        oppic_dat node_volume          = oppic_decl_dat(nodes_set, 1, OPP_TYPE_REAL, (char*)mesh.node_volume,  "node_volume");        
        oppic_dat node_potential       = oppic_decl_dat(nodes_set, 1, OPP_TYPE_REAL, (char*)mesh.node_pot,     "node_potential");     
        oppic_dat node_charge_density  = oppic_decl_dat(nodes_set, 1, OPP_TYPE_REAL, (char*)mesh.node_ion_den, "node_charge_density");
// opp_printf("XXXXXXXXXX", OPP_my_rank, " 65");
        // oppic_dat iface_v_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh.iface_v_normal,      "iface_v_normal");        
        // oppic_dat iface_u_normal       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh.iface_u_normal,      "iface_u_normal"); 
        // oppic_dat iface_normal         = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh.iface_normal,        "iface_normal");     
        // oppic_dat iface_area           = oppic_decl_dat(ifaces_set, 1,          OPP_TYPE_REAL, (char*)mesh.iface_area,          "iface_area");
        // oppic_dat iface_inj_part_dist  = oppic_decl_dat(ifaces_set, 1,          OPP_TYPE_INT,  (char*)mesh.iface_inj_part_dist, "iface_inj_part_dist");
        // oppic_dat iface_node_pos       = oppic_decl_dat(ifaces_set, DIMENSIONS, OPP_TYPE_REAL, (char*)mesh.iface_node_pos,      "iface_node_pos"); 
// opp_printf("XXXXXXXXXX", OPP_my_rank, " 66");
        oppic_set particles_set        = oppic_decl_particle_set(pcount, "particles", cells_set); 
        oppic_dat part_position        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     OPP_TYPE_REAL, (char*)arr_part_position, "part_position");
        oppic_dat part_velocity        = oppic_decl_particle_dat(particles_set, DIMENSIONS,     OPP_TYPE_REAL, (char*)arr_part_velocity, "part_velocity");    
        oppic_dat part_lc              = oppic_decl_particle_dat(particles_set, NODES_PER_CELL, OPP_TYPE_INT, (char*)arr_part_lc, "part_lc");
        oppic_dat part_mesh_relation   = oppic_decl_particle_dat(particles_set, 1,              OPP_TYPE_INT,  (char*)arr_part_mesh_relation, "part_mesh_relation", true); // new cell index field

delete[] arr_part_position;
delete[] arr_part_velocity;
delete[] arr_part_lc;
delete[] arr_part_mesh_relation;
    // TODO : Need to enrich mesh.dummy_part_random before using dummy_part
        // oppic_set dummy_part_set       = oppic_decl_particle_set(mesh.n_approx_injected, "dummy particles", cells_set); 
        // oppic_dat dum_part_random      = oppic_decl_dat(dummy_part_set, 2, OPP_TYPE_REAL, (char*)mesh.dummy_part_random, "dum_part_random");

// opp_printf("XXXXXXXXXX", OPP_my_rank, " 9");

        oppic_decl_const<OPP_REAL>(1, &spwt,         "CONST_spwt");
        oppic_decl_const<OPP_REAL>(1, &ion_velocity, "CONST_ion_velocity");
        oppic_decl_const<OPP_REAL>(1, &dt,           "CONST_dt");
        oppic_decl_const<OPP_REAL>(1, &plasma_den,   "CONST_plasma_den");
        oppic_decl_const<OPP_REAL>(1, &mass,         "CONST_mass");
        oppic_decl_const<OPP_REAL>(1, &charge,       "CONST_charge");

// opp_printf("XXXXXXXXXX", OPP_my_rank, " 8");

        mesh.DeleteValues();

// opp_printf("XXXXXXXXXX", OPP_my_rank, " 10");

        opp_partition(cells_set, cell_to_nodes_map);

//  opp_printf("XXXXXXXXXX", OPP_my_rank, " 67");

    if (false)
    {
        std::string f = std::string("A_") + std::to_string(OPP_my_rank);
        // oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
        // oppic_print_map_to_txtfile(cell_to_cell_map   , f.c_str(), "cell_to_cell_map.dat");
        // oppic_print_map_to_txtfile(iface_to_cell_map  , f.c_str(), "iface_to_cell_map.dat");

        // oppic_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");
        // oppic_print_dat_to_txtfile(node_potential     , f.c_str(), "node_potential.dat");
        // oppic_print_dat_to_txtfile(cell_electric_field, f.c_str(), "cell_electric_field.dat");

        // oppic_print_dat_to_txtfile(iface_v_normal,      f.c_str(), "iface_v_normal.dat");
        // oppic_print_dat_to_txtfile(iface_u_normal,      f.c_str(), "iface_u_normal.dat");
        // oppic_print_dat_to_txtfile(iface_normal,        f.c_str(), "iface_normal.dat");
        // oppic_print_dat_to_txtfile(iface_area,          f.c_str(), "iface_area.dat");
        // oppic_print_dat_to_txtfile(iface_inj_part_dist, f.c_str(), "iface_inj_part_dist.dat");

        oppic_print_dat_to_txtfile(part_position     , f.c_str(), "part_position.dat");
        oppic_print_dat_to_txtfile(part_velocity     , f.c_str(), "part_velocity.dat");
        oppic_print_dat_to_txtfile(part_lc           , f.c_str(), "part_lc.dat");
        oppic_print_dat_to_txtfile(part_mesh_relation, f.c_str(), "part_mesh_relation.dat");
    }

    if (true)
    {
        oppic_set set = particles_set;
        int comm_iteration = 0;
        int start1 = 0;
        int end = set->size;

        do
        {    
            opp_printf("Main()", "Starting iteration %d", comm_iteration);

            oppic_init_particle_move(set);

            int *mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

            for (int i = start1; i < end; i++)
            {        
                opp_move_var m;

                opp_printf("Main()", "Iter index %d comm_iteration %d *******************************", i, comm_iteration);

                do
                { 
                    int& map0idx      = mesh_relation_data[i];

                    // Expect this as the kernel **************************************************************
                    if (OPP_my_rank == 0 && i == 2 && comm_iteration==0)
                    {
                        opp_printf("EXPECT COMM", "3468 cell of INDEX 1 ieh of RANK 0");
                        ((int*)part_mesh_relation->data)[i] = 3862;
                        // ((int*)part_mesh_relation->data)[i] = 1800; 2457;
                        m.OPP_move_status = OPP_NEED_MOVE;
                    }
                    else if (OPP_my_rank == 0 && i == 19 && comm_iteration==0)
                    {
                        opp_printf("EXPECT COMM", "3649 cell of INDEX 182 ieh of RANK 0");
                        ((int*)part_mesh_relation->data)[i] = 4043;
                        m.OPP_move_status = OPP_NEED_MOVE;
                    }
                    else if (OPP_my_rank == 1 && i == 3 && comm_iteration==0)
                    {
                        opp_printf("EXPECT COMM", "3670 cell of INDEX 3 ieh of RANK 1");
                        ((int*)part_mesh_relation->data)[i] = 3653;
                        m.OPP_move_status = OPP_NEED_MOVE;
                    }
                    else if (OPP_my_rank == 1 && i == 19 && comm_iteration==1)
                    {
                        opp_printf("EXPECT COMM", "3672 cell of INDEX 1 ieh of RANK 0");
                        ((int*)part_mesh_relation->data)[i] = 3655;
                        m.OPP_move_status = OPP_NEED_MOVE;
                    }
                    else
                    {
                        m.OPP_move_status = OPP_MOVE_DONE;
                    }
                    // End of the kernel **************************************************************

                    // should check whether map0idx is in halo list, if yes, pack the particle data into MPI buffer
                    opp_check_part_need_comm(map0idx, set, i, m);

                } while (m.OPP_move_status == OPP_NEED_MOVE);

                if (m.OPP_move_status == OPP_NEED_REMOVE) 
                {
                    set->particle_remove_count += 1;
                    mesh_relation_data[i] = MAX_CELL_INDEX;
                }
            }

            if (oppic_finalize_particle_move(set))
            {
                // all mpi ranks do not have anything to communicate to any rank
                break;
            }
            else
            {
                // wait till all the particles are communicated and added to the dats
                opp_wait_all_particles(set);
                start1 = set->size - set->diff;
                end = set->size;
            }

            comm_iteration++;

        } while (true);
    }

if (true)
{
    std::string f = std::string("F_") + std::to_string(OPP_my_rank);

    oppic_print_dat_to_txtfile(part_position     , f.c_str(), "part_position.dat");
    oppic_print_dat_to_txtfile(part_velocity     , f.c_str(), "part_velocity.dat");
    oppic_print_dat_to_txtfile(part_lc           , f.c_str(), "part_lc.dat");
    oppic_print_dat_to_txtfile(part_mesh_relation, f.c_str(), "part_mesh_relation.dat");
}









/*


opp_printf("XXXXXXXXXX", OPP_my_rank, " 11");

        auto start = std::chrono::system_clock::now();
        auto start_iter1 = std::chrono::system_clock::now();
max_iter = 0;
        for (ts = 0; ts < max_iter; ts++)
        {
            if (ts == 1) start_iter1 = std::chrono::system_clock::now();

            int injected_count = 0;
        //     oppic_inject__Increase_particle_count(
        //         particles_set,                                                                // inject to particles_set
        //         ifaces_set,                                                                   // inlect_face_set
        //         oppic_arg_gbl(&(injected_count), 1, "int",   OP_RW),                          // injected total global,
        //         oppic_arg_dat(iface_area,                    OP_READ),                        // iface_area,
        //         oppic_arg_dat(iface_inj_part_dist,           OP_WRITE, OPP_Map_to_Mesh_Rel),  // iface_inj_part_dist,
        //         oppic_arg_gbl(&(remainder),     1, "double", OP_RW)                           // remainder global,
        //     );

        //     int old_nparts = particles_set->size;
        //     oppic_par_loop_inject__InjectIons(
        //         particles_set,                                                                               // particles_set
        //         oppic_arg_dat(part_position,                             OP_WRITE),                          // part_position,
        //         oppic_arg_dat(part_velocity,                             OP_WRITE),                          // part_velocity,
        //         oppic_arg_dat(part_mesh_relation,                        OP_RW),                             // part_cell_connectivity,
        //         oppic_arg_dat(iface_to_cell_map,                         OP_READ, OPP_Map_from_Mesh_Rel),    // iface to cell map
        //         oppic_arg_dat(cell_electric_field, 0, iface_to_cell_map, OP_READ, OPP_Map_from_Mesh_Rel),    // cell_ef,
        //         oppic_arg_dat(iface_u_normal,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_u,
        //         oppic_arg_dat(iface_v_normal,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_v,
        //         oppic_arg_dat(iface_normal,                              OP_READ, OPP_Map_from_Mesh_Rel),    // iface_normal,
        //         oppic_arg_dat(iface_node_pos,                            OP_READ, OPP_Map_from_Mesh_Rel),    // iface_node_pos
        //         oppic_arg_dat(dum_part_random,                           OP_READ, OPP_Map_from_Inj_part)     // dum_part_random
        //     );

            oppic_reset_dat(node_charge_density, (char*)opp_zero_double16);
        //     oppic_par_loop_particle_all__MoveToCells(
        //         particles_set,                                                                                // particles_set
        //         oppic_arg_dat(cell_electric_field,                       OP_READ, OPP_Map_from_Mesh_Rel),     // cell_ef,
        //         oppic_arg_dat(part_position,                             OP_RW),                              // part_pos,
        //         oppic_arg_dat(part_velocity,                             OP_RW),                              // part_vel,
        //         oppic_arg_dat(part_lc,                                   OP_RW),                              // part_lc,
        //         oppic_arg_dat(part_mesh_relation,                        OP_RW),                              // current_cell_index,
        //         oppic_arg_dat(cell_volume,                               OP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_volume,
        //         oppic_arg_dat(cell_determinants,                         OP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_det,
        //         oppic_arg_dat(cell_to_cell_map,                          OP_READ, OPP_Map_from_Mesh_Rel),     // cell_connectivity,
        //         oppic_arg_dat(node_charge_density, 0, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den0,
        //         oppic_arg_dat(node_charge_density, 1, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den1,
        //         oppic_arg_dat(node_charge_density, 2, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den2,
        //         oppic_arg_dat(node_charge_density, 3, cell_to_nodes_map, OP_INC,  OPP_Map_from_Mesh_Rel)      // node_charge_den3,
        //     );

            oppic_par_loop_all__ComputeNodeChargeDensity(
                nodes_set,                                       // nodes_set
                oppic_arg_dat(node_charge_density,  OP_RW),      // node_charge_density
                oppic_arg_dat(node_volume,          OP_READ)     // node_volume
            );

        //     mesh.solver->computePhi(  // TODO: Change this to kernel calls
        //         mesh.fesolver_method,
        //         oppic_arg_dat(node_potential,      OP_WRITE),
        //         oppic_arg_dat(node_charge_density, OP_READ)
        //     );

        //     oppic_reset_dat(cell_electric_field, (char*)opp_zero_double16); 
        //     oppic_par_loop_all__ComputeElectricField(
        //         cells_set,                                                        // cells_set
        //         oppic_arg_dat(cell_electric_field,                  OP_INC),      // cell_electric_field,
        //         oppic_arg_dat(cell_shape_deriv,                     OP_READ),     // cell_shape_deriv,
        //         oppic_arg_dat(node_potential, 0, cell_to_nodes_map, OP_READ),     // node_potential0,
        //         oppic_arg_dat(node_potential, 1, cell_to_nodes_map, OP_READ),     // node_potential1,
        //         oppic_arg_dat(node_potential, 2, cell_to_nodes_map, OP_READ),     // node_potential2,
        //         oppic_arg_dat(node_potential, 3, cell_to_nodes_map, OP_READ)      // node_potential3,
        //     );

        //     {
        //         double max_den = 0.0, max_phi = 0.0;
        //         if (params.get<OPP_BOOL>("check_max_values"))
        //         {
        //             for (int n = 0; n< mesh.n_nodes; n++) // ideally, need to copy data from device to host, but at this point host has correct data
        //             {
        //                 if (((double*)node_charge_density->data)[n] > max_den) max_den = ((double*)node_charge_density->data)[n];
        //                 if (abs(((double*)node_potential->data)[n]) > max_phi) max_phi = abs(((double*)node_potential->data)[n]);   
        //             }
        //         }

        //         std::cout<<"ts: " << ts << "\t np: " << particles_set->size 
        //             << " (" <<  injected_count << " added, "<< old_nparts - particles_set->size << " removed)"                
        //             << "\t max den: " << max_den << "\t max |phi|: " << max_phi << "\t ****" << std::endl;
        //     }
        }

        std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
        std::chrono::duration<double> diff_1 = std::chrono::system_clock::now() - start_iter1;
        opp_printf("FEMPIC", "Time to iterate ts=%d <sec>: ", max_iter);
        // std::cout << "\nFEMPIC - Time to iterate ts=" << max_iter << " <sec>: " << diff.count() << " and (ts-1) <sec>:" << diff_1.count() << "\n\n";
*/
        oppic_exit(); 
    // } // End Scope for oppic

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// oppic_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
// oppic_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");