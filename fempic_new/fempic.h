#pragma once

#include <oppic_lib.h>
#include <memory>
#include <regex>

#include "fempic_ori/meshes.h"
#include "fempic_ori/particles.h"
#include "fempic_ori/maths.h"
#include "FESolver.h"

using namespace std;

#define DIMENSIONS         3
#define NODES_PER_CELL     4
#define NEIGHBOUR_CELLS    4
#define DET_FIELDS         4
#define PRINT_PRECISION    15

#define MAX_CELL_INDEX     INT_MAX

class FieldPointers
{
    public:
        FieldPointers() {}
        virtual ~FieldPointers()
        {
            DeleteValues();   
        }

        inline void DeleteValues()
        {
            if (cell_to_nodes) delete[] cell_to_nodes;
            if (cell_to_cell) delete[] cell_to_cell;
            if (cell_ef) delete[] cell_ef;
            if (cell_det) delete[] cell_det;
            if (cell_volume) delete[] cell_volume;
            if (node_bnd_pot) delete[] node_bnd_pot;
            if (node_pot) delete[] node_pot;
            if (node_ion_den) delete[] node_ion_den;
            if (node_pos) delete[] node_pos;
            if (node_volume) delete[] node_volume;
            if (iface_to_cell) delete[] iface_to_cell;
            if (iface_to_nodes) delete[] iface_to_nodes;
            if (iface_v_normal) delete[] iface_v_normal;
            if (iface_u_normal) delete[] iface_u_normal;
            if (iface_normal) delete[] iface_normal;
            if (iface_area) delete[] iface_area;
            if (iface_inj_part_dist) delete[] iface_inj_part_dist;
            if (iface_node_pos) delete[] iface_node_pos;
            if (cell_shape_deriv) delete[] cell_shape_deriv;
            if (dummy_part_random) delete[] dummy_part_random;

            cell_to_nodes = nullptr;
            cell_to_cell = nullptr;
            cell_ef = nullptr;
            cell_det = nullptr;
            cell_volume = nullptr; 
            cell_shape_deriv = nullptr;

            node_bnd_pot = nullptr;
            node_pot = nullptr;
            node_ion_den = nullptr;
            node_pos = nullptr;
            node_volume = nullptr; 

            iface_to_cell = nullptr;    
            iface_to_nodes = nullptr;   
            iface_v_normal = nullptr;
            iface_u_normal = nullptr;
            iface_normal = nullptr;  
            iface_area = nullptr;    
            iface_inj_part_dist = nullptr;
            iface_node_pos = nullptr;
            dummy_part_random = nullptr;
        }

        int n_nodes;
        int n_cells;
        int n_ifaces;
        int n_approx_injected;

        int *cell_to_nodes = nullptr;
        int *cell_to_cell = nullptr;
        double *cell_ef = nullptr;
        double *cell_det = nullptr;
        double *cell_volume = nullptr; 
        double *cell_shape_deriv = nullptr;

        double *node_bnd_pot = nullptr;
        double *node_pot = nullptr;
        double *node_ion_den = nullptr;
        double *node_pos = nullptr;
        double *node_volume = nullptr; 

        int *iface_to_cell = nullptr;         // cell_con
        int *iface_to_nodes = nullptr;        // con[3]; 
        double *iface_v_normal = nullptr;     // v[3]
        double *iface_u_normal = nullptr;     // u[3]
        double *iface_normal = nullptr;       // normal[3]
        double *iface_area = nullptr;         // area
        int *iface_inj_part_dist = nullptr;
        double *iface_node_pos = nullptr;

        double * dummy_part_random = nullptr;

        std::shared_ptr<FESolver> solver;

        FESolver::Method fesolver_method = FESolver::Method::GaussSeidel;
};

void oppic_seq_loop_inject__Increase_particle_count
(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
);

void oppic_par_loop_inject__InjectIons(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_velocity,
    oppic_arg arg2,     // part_cell_connectivity,
    oppic_arg arg3,     // iface to cell map
    oppic_arg arg4,     // cell_ef,
    oppic_arg arg5,     // iface_u,
    oppic_arg arg6,     // iface_v,
    oppic_arg arg7,     // iface_normal,
    oppic_arg arg8,     // iface_node_pos
    oppic_arg arg9      // dummy_part_random
);

void oppic_par_loop_particle_all__MoveToCells(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // cell_ef,
    oppic_arg arg1,     // part_pos,
    oppic_arg arg2,     // part_vel,
    oppic_arg arg3,     // part_lc,
    oppic_arg arg4,     // current_cell_index,
    oppic_arg arg5,     // current_cell_volume,
    oppic_arg arg6,     // current_cell_det,
    oppic_arg arg7,     // cell_connectivity,
    oppic_arg arg8,     // node_charge_den0,
    oppic_arg arg9,     // node_charge_den1,
    oppic_arg arg10,    // node_charge_den2,
    oppic_arg arg11     // node_charge_den3,
);

void oppic_par_loop_all__ComputeNodeChargeDensity(
    oppic_set set,     // nodes_set
    oppic_arg arg0,    // node_charge_density
    oppic_arg arg1     // node_volume
);

void oppic_par_loop_all__ComputeElectricField(
    oppic_set set,      // cells_set
    oppic_arg arg0,     // cell_electric_field,
    oppic_arg arg1,     // cell_shape_deriv,
    oppic_arg arg2,     // node_potential0,
    oppic_arg arg3,     // node_potential1,
    oppic_arg arg4,     // node_potential2,
    oppic_arg arg5      // node_potential3,
);

inline FieldPointers LoadMesh(opp::Params& params, int argc, char **argv)
{ TRACE_ME;

    FieldPointers mesh;

    std::shared_ptr<Volume> volume(new Volume());
    if (!LoadVolumeMesh(params.get<STRING>("global_mesh"), *(volume.get())) ||
        !LoadSurfaceMesh(params.get<STRING>("inlet_mesh"), *(volume.get()),INLET, params.get<BOOL>("invert_normals")) ||
        !LoadSurfaceMesh(params.get<STRING>("wall_mesh"), *(volume.get()), FIXED, params.get<BOOL>("invert_normals"))) return mesh;

    volume->summarize(std::cout);

    mesh.n_nodes         = volume->nodes.size();
    mesh.n_cells         = volume->elements.size();
    mesh.n_ifaces        = volume->inlet_faces.size();

    mesh.cell_ef         = new double[mesh.n_cells * DIMENSIONS];
    mesh.cell_to_nodes   = new int[mesh.n_cells * NODES_PER_CELL];
    mesh.cell_to_cell    = new int[mesh.n_cells * NEIGHBOUR_CELLS];
    mesh.cell_det        = new double[mesh.n_cells * DET_FIELDS * NEIGHBOUR_CELLS]; // arranged as [alpha,beta,gamma,delta] * 4 neighbours
    mesh.cell_volume     = new double[mesh.n_cells];
    mesh.cell_shape_deriv = new double[mesh.n_cells * NODES_PER_CELL*DIMENSIONS]; // arranged as [x,y,z] * 4 neighbours

    mesh.node_bnd_pot    = new double[mesh.n_nodes];
    mesh.node_pot        = new double[mesh.n_nodes];
    mesh.node_ion_den    = new double[mesh.n_nodes];
    mesh.node_pos        = new double[mesh.n_nodes * DIMENSIONS];
    mesh.node_volume     = new double[mesh.n_nodes];

    mesh.iface_to_cell   = new int[mesh.n_ifaces];
    mesh.iface_to_nodes  = new int[mesh.n_ifaces * DIMENSIONS]; 
    mesh.iface_v_normal  = new double[mesh.n_ifaces * DIMENSIONS]; 
    mesh.iface_u_normal  = new double[mesh.n_ifaces * DIMENSIONS]; 
    mesh.iface_normal    = new double[mesh.n_ifaces * DIMENSIONS]; 
    mesh.iface_area      = new double[mesh.n_ifaces]; 
    mesh.iface_inj_part_dist = new int[mesh.n_ifaces];
    mesh.iface_node_pos  = new double[mesh.n_ifaces * DIMENSIONS]; 

    for (int n=0; n<mesh.n_nodes; n++)
    {
        switch (volume->nodes[n].type)
        {
            case INLET:    mesh.node_bnd_pot[n] = 0; break;                                         /*phi_inlet*/
            case FIXED:    mesh.node_bnd_pot[n] = -(params.get<REAL>("wall_potential")); break;      /*fixed phi points*/
            default:       mesh.node_bnd_pot[n] = 0;                                                /*default*/
        }

        mesh.node_pot[n] = 0.0f;
        mesh.node_ion_den[n] = 0.0f;

        Node &node = volume->nodes[n];
    
        for (int dim=0; dim<DIMENSIONS; dim++)
            mesh.node_pos[n * DIMENSIONS + dim] = node.pos[dim];
    
        mesh.node_volume[n]     = node.volume;
    }

    for (int cellID=0; cellID<mesh.n_cells; cellID++)
    {
        Tetra &tet = volume->elements[cellID];
        
        for (int nodeCon=0; nodeCon<NODES_PER_CELL; nodeCon++)
        {
            mesh.cell_to_nodes[cellID * NODES_PER_CELL + nodeCon] = tet.con[nodeCon];

            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 0 ] = 0.0;
            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 1 ] = 0.0;
            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 2 ] = 0.0;
        }
        
        for (int cellCon=0; cellCon<NEIGHBOUR_CELLS; cellCon++)
        {
            mesh.cell_to_cell[cellID * NEIGHBOUR_CELLS + cellCon]     = tet.cell_con[cellCon];

            mesh.cell_det[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 0] = tet.alpha[cellCon];
            mesh.cell_det[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 1] = tet.beta[cellCon];
            mesh.cell_det[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 2] = tet.gamma[cellCon];
            mesh.cell_det[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 3] = tet.delta[cellCon];
        }

        mesh.cell_volume[cellID] = tet.volume;

        mesh.cell_ef[cellID * DIMENSIONS + 0] = 0.0;
        mesh.cell_ef[cellID * DIMENSIONS + 1] = 0.0;
        mesh.cell_ef[cellID * DIMENSIONS + 2] = 0.0;
    }

    for (int faceID=0; faceID<mesh.n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        Node &node = volume->nodes[face.con[0]];

        mesh.iface_to_cell[faceID] = face.cell_con;
        mesh.iface_area[faceID] = face.area;
        mesh.iface_inj_part_dist[faceID] = 0;

        for (int dim=0; dim<3; dim++)
        {
            mesh.iface_to_nodes[faceID * 3 + dim] = face.con[dim]; 
            mesh.iface_v_normal[faceID * 3 + dim] = face.v[dim];   
            mesh.iface_u_normal[faceID * 3 + dim] = face.u[dim];   
            mesh.iface_normal[faceID * 3 + dim]   = face.normal[dim]; 

            mesh.iface_node_pos[faceID * 3 + dim] = node.pos[dim];
        }
    }

    mesh.solver = std::make_shared<FESolver>(volume, argc, argv);

    mesh.solver->phi0 = 0;
    mesh.solver->n0 = params.get<REAL>("plasma_den");
    mesh.solver->kTe = Kb * params.get<REAL>("electron_temperature");

    for (int n = 0; n < mesh.n_nodes; n++) // This g array is a duplicate, can remove and attach node_bnd_pot dat instead
    {
        switch (volume->nodes[n].type)
        {
            case INLET:    mesh.solver->g[n] = 0; break;                                         /*phi_inlet*/
            case FIXED:    mesh.solver->g[n] = -(params.get<REAL>("wall_potential")); break;     /*fixed phi points*/
            default:       mesh.solver->g[n] = 0;                                                /*default*/
        }
    }

    mesh.solver->startAssembly();
    mesh.solver->preAssembly();    /*this will form K and F0*/

    mesh.solver->summarize(std::cout);

    if      (std::regex_match(params.get<STRING>("fesolver_method"), std::regex("nonlinear", std::regex_constants::icase))) mesh.fesolver_method = FESolver::NonLinear;
    else if (std::regex_match(params.get<STRING>("fesolver_method"), std::regex("gaussseidel", std::regex_constants::icase))) mesh.fesolver_method = FESolver::GaussSeidel;
    else if (std::regex_match(params.get<STRING>("fesolver_method"), std::regex("lapack", std::regex_constants::icase))) mesh.fesolver_method = FESolver::Lapack;
    else if (std::regex_match(params.get<STRING>("fesolver_method"), std::regex("petsc", std::regex_constants::icase))) mesh.fesolver_method = FESolver::Petsc;

    for (int cellID=0; cellID<mesh.n_cells; cellID++)
    {
        for (int nodeCon=0; nodeCon<NODES_PER_CELL; nodeCon++)
        {
            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 0 ] = mesh.solver->NX[cellID][nodeCon][0];
            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 1 ] = mesh.solver->NX[cellID][nodeCon][1];
            mesh.cell_shape_deriv[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 2 ] = mesh.solver->NX[cellID][nodeCon][2];
        }
    }

    double plasma_den = params.get<REAL>("plasma_den");
    double dt = params.get<REAL>("dt");
    double ion_velocity = params.get<REAL>("ion_velocity");
    double spwt = params.get<REAL>("spwt");

    mesh.n_approx_injected = 0;
    double remainder = 0.0;

    for (int faceID=0; faceID<mesh.n_ifaces; faceID++)
    {   
        {   // DUPLICATE: This calculation is in kernels
            double num_per_sec = plasma_den * ion_velocity * mesh.iface_area[faceID];
            double num_real = num_per_sec * dt;
            double fnum_mp = num_real / spwt + remainder;
            int num_mp = (int)fnum_mp;
            remainder = fnum_mp - num_mp;
            mesh.n_approx_injected += num_mp;
        }
    }

    mesh.n_approx_injected += 100; // Add a few (100) for safety
    
    mesh.dummy_part_random  = new double[mesh.n_approx_injected * 2]; 

    for (int i=0; i<mesh.n_approx_injected; i++)
    {   
        mesh.dummy_part_random[i * 2 + 0] = rnd();
        mesh.dummy_part_random[i * 2 + 1] = rnd();
    }

    return mesh;
}
