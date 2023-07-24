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

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

// #ifndef ENABLE_MPI
//     #define ENABLE_MPI
// #endif

#include "opp_move_with_approx.h"
extern std::unique_ptr<opp::CellApproximator> opp_mover;

#include <oppic_lib.h>
#ifdef ENABLE_MPI
    #include <opp_mpi.h>
#endif

#include "fempic_ori/meshes.h"
#include "fempic_ori/particles.h"
#include "fempic_ori/maths.h"
#include "FESolver.h"
#include "cluster.h"

#define USE_RAND_FILE
#define USE_CELL_PARTITIONING

#ifdef DEBUG_LOG
    #define FP_DEBUG true
#else
    #define FP_DEBUG false
#endif

#define ONE                1
#define DIM                3
#define N_PER_C            4
#define NEIGHB_C           4
#define DET_FIELDS         4
#define ALL_DET            (NEIGHB_C * DET_FIELDS)
#define PRINT_PRECISION    15

#define MAX_c_INDEX        INT_MAX
#define INJ_EXCESS         0 // 100

// #define DET_SCALE          1.0
// #define SCALE2             1.0

// extern int ts;
// using namespace std;

//*************************************************************************************************
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
            if (c_to_n) delete[] c_to_n;
            if (c_to_c) delete[] c_to_c;
            if (c_ef) delete[] c_ef;
            if (c_det) delete[] c_det;
            if (c_vol) delete[] c_vol;
            if (c_sd) delete[] c_sd;
            if (c_col) delete[] c_col;
            if (c_id) delete[] c_id;
            if (n_bnd_pot) delete[] n_bnd_pot;
            if (n_pot) delete[] n_pot;
            if (n_ion_den) delete[] n_ion_den;
            if (n_pos) delete[] n_pos;
            if (n_vol) delete[] n_vol;
            if (n_type) delete[] n_type;
            if (n_color) delete[] n_color;
            if (n_id) delete[] n_id;
            if (if_to_c) delete[] if_to_c;
            if (if_to_n) delete[] if_to_n;
            if (if_v_norm) delete[] if_v_norm;
            if (if_u_norm) delete[] if_u_norm;
            if (if_norm) delete[] if_norm;
            if (if_area) delete[] if_area;
            if (if_dist) delete[] if_dist;
            if (if_n_pos) delete[] if_n_pos;
            if (if_id) delete[] if_id;

            c_to_n = nullptr;
            c_to_c = nullptr;
            c_ef = nullptr;
            c_det = nullptr;
            c_vol = nullptr; 
            c_sd = nullptr;
            c_col = nullptr;
            c_id = nullptr;

            n_bnd_pot = nullptr;
            n_pot = nullptr;
            n_ion_den = nullptr;
            n_pos = nullptr;
            n_vol = nullptr; 
            n_type = nullptr;
            n_color = nullptr;
            n_id = nullptr;

            if_to_c = nullptr;    
            if_to_n = nullptr;   
            if_v_norm = nullptr;
            if_u_norm = nullptr;
            if_norm = nullptr;  
            if_area = nullptr;    
            if_dist = nullptr;
            if_n_pos = nullptr;
            if_id = nullptr;
        }

        inline void CreateMeshArrays()
        {
            c_ef      = new double[n_cells * DIM];
            c_to_n    = new int[n_cells * N_PER_C];
            c_to_c    = new int[n_cells * NEIGHB_C];
            c_det     = new double[n_cells * ALL_DET]; // [alpha,beta,gamma,delta] * 4neighbours
            c_vol     = new double[n_cells];
            c_sd      = new double[n_cells * N_PER_C*DIM]; // arranged as [x,y,z] * 4 neighbours
            c_col     = new int[n_cells];
            c_id      = new int[n_cells];

            n_bnd_pot = new double[n_nodes];
            n_pot     = new double[n_nodes];
            n_ion_den = new double[n_nodes];
            n_pos     = new double[n_nodes * DIM];
            n_vol     = new double[n_nodes];
            n_type    = new int[n_nodes];
            n_color   = new int[n_nodes];
            n_id      = new int[n_nodes];

            if_to_c   = new int[n_ifaces];
            if_to_n   = new int[n_ifaces * DIM]; 
            if_v_norm = new double[n_ifaces * DIM]; 
            if_u_norm = new double[n_ifaces * DIM]; 
            if_norm   = new double[n_ifaces * DIM]; 
            if_area   = new double[n_ifaces]; 
            if_dist   = new int[n_ifaces];
            if_n_pos  = new double[n_ifaces * DIM];
            if_id     = new int[n_ifaces];  
        };

        int n_nodes  = 0;
        int n_cells  = 0;
        int n_ifaces = 0;

        int *c_to_n = nullptr;
        int *c_to_c = nullptr;
        double *c_ef = nullptr;
        double *c_det = nullptr;
        double *c_vol = nullptr; 
        double *c_sd = nullptr;
        int *c_col = nullptr; 
        int *c_id = nullptr;

        double *n_bnd_pot = nullptr;
        double *n_pot = nullptr;
        double *n_ion_den = nullptr;
        double *n_pos = nullptr;
        double *n_vol = nullptr; 
        int *n_type = nullptr; 
        int *n_color = nullptr; 
        int *n_id = nullptr;

        int *if_to_c = nullptr;        // c_con
        int *if_to_n = nullptr;        // con[3]; 
        double *if_v_norm = nullptr;   // v[3]
        double *if_u_norm = nullptr;   // u[3]
        double *if_norm = nullptr;     // normal[3]
        double *if_area = nullptr;     // area
        int *if_dist = nullptr;
        double *if_n_pos = nullptr;
        int *if_id = nullptr;
};

//*************************************************************************************************
inline int InitializeInjectDistributions(oppic_dat if_dist_dat, oppic_dat if_area_dat, 
                                        oppic_dat dummy_random)
{
    if (OP_DEBUG) opp_printf("InitializeInjectDistributions", "START");

    double plasma_den   = opp_params->get<OPP_REAL>("plasma_den");
    double dt           = opp_params->get<OPP_REAL>("dt");
    double ion_velocity = opp_params->get<OPP_REAL>("ion_velocity");
    double spwt         = opp_params->get<OPP_REAL>("spwt");

    int total_inject_count = 0, max_inject_count_per_face = 0;
    double remainder = 0.0;

    // find the number of particles to be injected through each inlet face and 
    // get the max injected particle count per face
    for (int faceID=0; faceID<if_area_dat->set->size; faceID++)
    {   
        remainder = 0.0; // always make remainder to zero, and if not the MPI results will change

        double num_per_sec = plasma_den * ion_velocity * ((double*)if_area_dat->data)[faceID];
        double num_real = num_per_sec * dt;
        double fnum_mp = num_real / spwt + remainder;
        int num_mp = (int)fnum_mp;
        remainder = fnum_mp - num_mp;

        total_inject_count += num_mp;

        ((int*)if_dist_dat->data)[faceID] = total_inject_count;

        if (max_inject_count_per_face < num_mp)
            max_inject_count_per_face = num_mp;
    }

    // increase dummy random particle set size, to load the random numbers for particle injections
    oppic_increase_particle_count(dummy_random->set, (total_inject_count + INJ_EXCESS));

#ifdef USE_RAND_FILE
    
    if (OP_DEBUG)
        opp_printf("InitializeInjectDistributions", "RAND_FILE inj_count %d max_inj_count_per_face %d", 
            total_inject_count, max_inject_count_per_face);    

    int total_size = -1, fsize = -1, fdim = -1;
    FILE *fp = NULL;
    std::string rand_file_path = opp_params->get<OPP_STRING>("rand_file");

    // read from MPI ROOT and broadcast // alternatively could use HDF5 files
    if (OPP_rank == OPP_ROOT)
    {       
        if ((fp = fopen(rand_file_path.c_str(), "r")) == NULL)
        {
            opp_printf("InitializeInjectDistributions", "Unable to open file %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2)
        {
            opp_printf("InitializeInjectDistributions", "Error reading file data from %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        total_size = fsize * fdim;

        if (max_inject_count_per_face * dummy_random->dim > total_size)
        {
            opp_printf("InitializeInjectDistributions", "dim and/or set_size issue in file %s\n", 
                rand_file_path.c_str());
            opp_abort();     
        }
    }

#ifdef ENABLE_MPI
    // Load the whole file and bradcast
    // TODO : We can reduce communications by sending only the required size
    MPI_Bcast(&total_size, 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
#endif

    double* dist = new double[total_size];
    if (OPP_rank == OPP_ROOT)
    {
        for (int n = 0; n < fsize; n++)
        {
            if (fscanf(fp, " %lf %lf\n", &dist[n * 2 + 0], &dist[n * 2 + 1]) != 2) 
            {
                opp_printf("InitializeInjectDistributions", "Error reading from %s at index %d\n", 
                    rand_file_path.c_str(), n);
                opp_abort();
            }
        }

        fclose(fp);
    }

#ifdef ENABLE_MPI
    // Load the whole file and bradcast
    // TODO : We can reduce communications by sending only the required size
    MPI_Bcast(dist, total_size, MPI_DOUBLE, OPP_ROOT, MPI_COMM_WORLD);
#endif

    double* random_dat = (double *)dummy_random->data;
    int* distribution  = (int *)if_dist_dat->data;
    int iface_index = 0, part_in_face_index = 0;

    // This trouble is only because we need mpi results to match with seq and others
    for (int i = 0; i < dummy_random->set->size; i++)
    {
        if (i >= distribution[iface_index])
        {
            iface_index++; // check whether it is j or j-1
            part_in_face_index = 0; 
        }

        random_dat[i * 2 + 0] = dist[part_in_face_index * 2 + 0];
        random_dat[i * 2 + 1] = dist[part_in_face_index * 2 + 1];

        part_in_face_index++;
    }

    delete[] dist;

#else
    // Might not be able to verify MPI results with seq with this

    if (OP_DEBUG)
        opp_printf("InitializeInjectDistributions", "RAND_GENERATE inj_count %d max_inj_count_per_face %d", 
            total_inject_count, max_inject_count_per_face);    
    
    double *dist = (double*)dummy_random->data;

    for (int i = 0; i < dummy_random->set->size * dummy_random->dim; i++)
    {   
        dist[i] = rnd();
    }   

#endif

    dummy_random->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    if_dist_dat->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    
    opp_printf("InitializeInjectDistributions", "total_inject_count %d", total_inject_count);

    return total_inject_count;
}

//*************************************************************************************************
inline std::shared_ptr<FieldPointers> LoadMesh()
{
    std::shared_ptr<FieldPointers> mesh(new FieldPointers());

    std::shared_ptr<Volume> volume(new Volume());
    if (!LoadVolumeMesh(opp_params->get<OPP_STRING>("global_mesh"), *(volume.get())) ||
        !(LoadSurfaceMesh(opp_params->get<OPP_STRING>("inlet_mesh"), *(volume.get()),INLET, 
            opp_params->get<OPP_BOOL>("invert_normals"))) ||
        !(LoadSurfaceMesh(opp_params->get<OPP_STRING>("wall_mesh"), *(volume.get()), FIXED, 
            opp_params->get<OPP_BOOL>("invert_normals")))) 
    {    
        return mesh;
    }

    if (OPP_rank == OPP_ROOT)
        volume->summarize(std::cout);

    mesh->n_nodes         = volume->nodes.size();
    mesh->n_cells         = volume->elements.size();
    mesh->n_ifaces        = volume->inlet_faces.size();

    mesh->CreateMeshArrays();

    for (int n=0; n<mesh->n_nodes; n++)
    {
        switch (volume->nodes[n].type)
        {
            case INLET: mesh->n_bnd_pot[n] = 0; break;
            case FIXED: mesh->n_bnd_pot[n] = -(opp_params->get<OPP_REAL>("wall_potential")); break;
            default:    mesh->n_bnd_pot[n] = 0; /*default*/
        }

        mesh->n_pot[n] = 0.0f;
        mesh->n_ion_den[n] = 0.0f;

        Node &node = volume->nodes[n];
    
        for (int dim=0; dim<DIM; dim++)
            mesh->n_pos[n * DIM + dim] = node.pos[dim];
    
        mesh->n_vol[n] = node.volume;
        mesh->n_type[n]   = (int)node.type;
        mesh->n_id[n]     = n;
    }

    for (int cID=0; cID<mesh->n_cells; cID++)
    {
        Tetra &tet = volume->elements[cID];
        
        for (int nodeCon=0; nodeCon<N_PER_C; nodeCon++)
        {
            mesh->c_to_n[cID * N_PER_C + nodeCon] = tet.con[nodeCon];

            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 0 ] = 0.0;
            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 1 ] = 0.0;
            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 2 ] = 0.0;
        }

        for (int cConn=0; cConn<NEIGHB_C; cConn++)
        {
            mesh->c_to_c[cID * NEIGHB_C + cConn]     = tet.cell_con[cConn];

            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 0] = (tet.alpha[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 1] = (tet.beta[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 2] = (tet.gamma[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 3] = (tet.delta[cConn]);
        }

        mesh->c_vol[cID] = tet.volume;
        mesh->c_id[cID] = cID;

        mesh->c_ef[cID * DIM + 0] = 0.0;
        mesh->c_ef[cID * DIM + 1] = 0.0;
        mesh->c_ef[cID * DIM + 2] = 0.0;        
    }

    for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        Node &node = volume->nodes[face.con[0]];

        mesh->if_to_c[faceID] = face.cell_con;
        mesh->if_area[faceID] = face.area;
        mesh->if_dist[faceID] = 0;
        mesh->if_id[faceID] = faceID;

        for (int dim=0; dim<3; dim++)
        {
            mesh->if_to_n[faceID * 3 + dim] = face.con[dim]; 
            mesh->if_v_norm[faceID * 3 + dim] = face.v[dim];   
            mesh->if_u_norm[faceID * 3 + dim] = face.u[dim];   
            mesh->if_norm[faceID * 3 + dim]   = face.normal[dim]; 

            mesh->if_n_pos[faceID * 3 + dim] = node.pos[dim];
        }
    }

#if defined(USE_NODE_PARTITIONING)

    // Cluster iface centroids and assign nodes to the rank on the major particle movement axis z
    std::vector<Point3D> face_points;
    
    for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        
        const Point3D *node0_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[0]].pos);
        const Point3D *node1_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[1]].pos);
        const Point3D *node2_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[2]].pos);

        face_points.push_back(getTriangleCentroid3D(*node0_pos, *node1_pos, *node2_pos));
    }

    std::vector<int> cluster_assignments = kMeansClustering3D(face_points, OPP_comm_size);

    std::vector<Point3D> cluster_centroids = calculateTriangleCentroids3D(face_points, cluster_assignments);

    // Print cluster assignments
    if (OP_DEBUG)
    {
        for (int i = 0; i < face_points.size(); ++i) {
            std::cout << "[" << face_points[i].x << "," << face_points[i].y << "," << 
                face_points[i].z << "], ";
        }
    
        std::cout << std::endl << std::endl;

        for (int i = 0; i < cluster_assignments.size(); ++i) {
            std::cout << cluster_assignments[i] << ",";
        }

        std::cout << std::endl << std::endl;

        for (int i = 0; i < cluster_centroids.size(); ++i) {
            std::cout << "[" << cluster_centroids[i].x << ","<< cluster_centroids[i].y << "," << 
                cluster_centroids[i].z << "], ";
        }
    }

    std::cout << std::endl << std::endl;
    bool color_found = false;

    for (int n=0; n<mesh->n_nodes; n++) {
        Node &node = volume->nodes[n];
        const Point3D *n_pos = reinterpret_cast<Point3D*>(node.pos);
        color_found = false;

        for (int faceID=0; faceID<mesh->n_ifaces; faceID++) {
            Face &face = volume->inlet_faces[faceID];
            
            const Point3D *node0_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[0]].pos);
            const Point3D *node1_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[1]].pos);
            const Point3D *node2_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[2]].pos);

            bool isInTriangle = isPointInTriangle(*n_pos, *node0_pos, *node1_pos, *node2_pos);
    
            if (isInTriangle) {
                mesh->n_color[n] = cluster_assignments[faceID];
                color_found = true;
                break;
            }
        }

        if (!color_found) {
            // if not found, then assign to the closest cluster centroid

            double min_idx = std::numeric_limits<double>::max();
            for (int i = 0; i < cluster_centroids.size(); ++i) {
                double d = euclideanDistancePlaneXY(cluster_centroids[i], *n_pos);
                if (d < min_idx) {
                    min_idx = d;
                    mesh->n_color[n] = cluster_assignments[i];
                }
            }

            // opp_printf("Setup", "Error... Couldnt find colour for node [%lf,%lf]", 
            //    n_pos->x, n_pos->y);
        }
    }

    int* counts = new int[OPP_comm_size];
    for (int n=0; n<OPP_comm_size; n++) 
        counts[n] = 0;

    for (int n=0; n<mesh->n_nodes; n++) {
        counts[mesh->n_color[n]]++;
    }

    if (OP_DEBUG)
    {
        for (int n=0; n<OPP_comm_size; n++) {
            opp_printf("Setup", "Rank|nodes %d %d", n, counts[n]);
        }
    }   

    Cluster iface centroids and assign cells to the rank on the major particle movement axis z

#elif defined(USE_CELL_PARTITIONING)

    std::vector<Point3D> face_points;
    
    for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        
        const Point3D *node0_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[0]].pos);
        const Point3D *node1_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[1]].pos);
        const Point3D *node2_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[2]].pos);

        face_points.push_back(getTriangleCentroid3D(*node0_pos, *node1_pos, *node2_pos));
    }

    std::vector<int> cluster_assignments = kMeansClustering3D(face_points, OPP_comm_size);

    std::vector<Point3D> cluster_centroids = calculateTriangleCentroids3D(face_points, cluster_assignments);

    bool color_found = false;
    std::vector<double> rank_volume(OPP_comm_size, 0.0);

    for (int cellID=0; cellID<mesh->n_cells; cellID++)
    {
        Tetra &tet = volume->elements[cellID];
        
        const Point3D *node0_pos = reinterpret_cast<Point3D*>(volume->nodes[tet.con[0]].pos);
        const Point3D *node1_pos = reinterpret_cast<Point3D*>(volume->nodes[tet.con[1]].pos);
        const Point3D *node2_pos = reinterpret_cast<Point3D*>(volume->nodes[tet.con[2]].pos);
        const Point3D *node3_pos = reinterpret_cast<Point3D*>(volume->nodes[tet.con[3]].pos);

        const Point3D c_pos = getTetraCentroid3D(*node0_pos, *node1_pos, *node2_pos, *node3_pos);

        color_found = false;

        for (int faceID=0; faceID<mesh->n_ifaces; faceID++) {
            Face &face = volume->inlet_faces[faceID];
            
            const Point3D *face_node0_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[0]].pos);
            const Point3D *face_node1_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[1]].pos);
            const Point3D *face_node2_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[2]].pos);

            bool isInTriangle = isPointInTriangle(c_pos, *face_node0_pos, *face_node1_pos, *face_node2_pos);

            if (isInTriangle) {
                mesh->c_col[cellID] = cluster_assignments[faceID];    
                color_found = true;
                rank_volume[cluster_assignments[faceID]] += (calculateTetraVolume(*node0_pos, *node1_pos, 
                    *node2_pos, *node3_pos) * 1000000000);
                break;
            }
        }

        if (!color_found) {
            opp_printf("Setup", "Error... Couldnt find colour for cell ,[%lf,%lf]", c_pos.x, c_pos.y);
        }   
    }

    // Print cluster assignments
    if (OP_DEBUG)
    {
        for (int i = 0; i < (int)face_points.size(); ++i) {
            std::cout << "[" << face_points[i].x << "," << face_points[i].y << "," << 
                face_points[i].z << "], ";
        }
    
        std::cout << std::endl << std::endl;

        for (int i = 0; i < (int)cluster_assignments.size(); ++i) {
            std::cout << cluster_assignments[i] << ",";
        }

        std::cout << std::endl << std::endl;

        for (int i = 0; i < (int)cluster_centroids.size(); ++i) {
            std::cout << "[" << cluster_centroids[i].x << ","<< cluster_centroids[i].y << "," << 
                cluster_centroids[i].z << "], ";
        }

        std::cout << std::endl << std::endl;

        // DEBUG ONLY - to be removed start
            std::vector<double> face_areas_per_rank(OPP_comm_size, 0.0);
            for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
            {
                Face &face = volume->inlet_faces[faceID];
                
                const Point3D *node0_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[0]].pos);
                const Point3D *node1_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[1]].pos);
                const Point3D *node2_pos = reinterpret_cast<Point3D*>(volume->nodes[face.con[2]].pos);

                face_areas_per_rank[cluster_assignments[faceID]] += 
                    std::abs(calculateTriangleArea(*node0_pos, *node1_pos, *node2_pos) * 1000000);
            }
            double total = 0.0;
            std::cout << "face_areas_per_rank (um^2) -> ";
            for (int i = 0; i < OPP_comm_size; ++i) {
                total += face_areas_per_rank[i];
                std::cout << i << "|" << face_areas_per_rank[i] << "|" << rank_volume[i] << " ";
            }
            std::cout << " = all|" << total << std::endl << std::endl;
        // DEBUG ONLY - to be removed end

        std::vector<int> c_counts(OPP_comm_size, 0);
        std::vector<int> face_counts(OPP_comm_size, 0);

        for (int n=0; n<mesh->n_cells; n++)
            c_counts[mesh->c_col[n]]++;

        for (int n=0; n<(int)cluster_assignments.size(); n++)
            face_counts[cluster_assignments[n]]++;

        std::string log = "";
        for (int n=0; n<OPP_comm_size; n++)
            log += std::to_string(n) + "|" + std::to_string(face_counts[n]) + "|" + \
                std::to_string(c_counts[n]) + " ";

        opp_printf("Setup", "Rank|faces|cells %s", log.c_str());
    }   

#endif

    return mesh;
}

//*************************************************************************************************
inline void print_per_cell_particle_counts(oppic_dat c_part_count, oppic_dat part_mesh_relation)     
{
    // validate the final particle counts in each cell
    opp_reset_dat(c_part_count, (char*)opp_zero_int16);
    for (int p = 0; p < part_mesh_relation->set->size; p++)
    {
        int c_index    = ((int *)part_mesh_relation->data)[p];
        ((int *)c_part_count->data)[c_index] += 1;
    }

#ifdef ENABLE_MPI
    opp_mpi_print_dat_to_txtfile(c_part_count, "c_part_count.dat");
#else
    opp_print_dat_to_txtfile(c_part_count, "", "c_part_count.dat");
#endif
}

//*************************************************************************************************
inline std::string get_global_level_log(oppic_dat n_charge_density, oppic_dat n_potential, int local_part_count, 
    int local_parts_injected, int local_part_removed)
{
    std::string log = "";
    double max_den = 0.0, max_phi = 0.0;
    double global_max_den = 0.0, global_max_phi = 0.0;
    int global_part_size = 0, global_inj_size = 0, global_removed = 0, global_max_comm_iteration = 0;

    // ideally, need to copy data from device to host, but at this point host has correct data
    for (int n = 0; n< n_potential->set->size; n++) 
    {
        if (abs(((double*)n_charge_density->data)[n]) > max_den) max_den = abs(((double*)n_charge_density->data)[n]);
        if (abs(((double*)n_potential->data)[n]) > max_phi) max_phi = abs(((double*)n_potential->data)[n]);   
    }

#ifdef ENABLE_MPI
    MPI_Reduce(&max_den, &global_max_den, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&max_phi, &global_max_phi, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&local_part_count, &global_part_size, 1, MPI_INT, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&local_parts_injected, &global_inj_size, 1, MPI_INT, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&local_part_removed, &global_removed, 1, MPI_INT, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&OPP_max_comm_iteration, &global_max_comm_iteration, 1, MPI_INT, MPI_MAX, OPP_ROOT, MPI_COMM_WORLD);
#else
    global_max_den = max_den;
    global_max_phi = max_phi;
    global_part_size = local_part_count;
    global_inj_size = local_parts_injected;
    global_removed = local_part_removed;
    global_max_comm_iteration = OPP_max_comm_iteration;
#endif

    log += std::string("\t np: ") + str(global_part_size, "%d");
    log += std::string(" (") + str(global_inj_size, "%d");
    log += std::string(" added, ") + str(global_removed, "%d");
    log += std::string(" removed)\t max den: ") + str(global_max_den, "%2.25lE");
    log += std::string(" max |phi|: ") + str(global_max_phi, "%2.10lE");
    log += std::string(" max_comm_iteration: ") + str(global_max_comm_iteration, "%d");
    return log;
}

//*************************************************************************************************
inline void DistributeMeshOverRanks(std::shared_ptr<FieldPointers>& g_m, std::shared_ptr<FieldPointers>& m)
{ 
#ifdef ENABLE_MPI
    MPI_Bcast(&(g_m->n_nodes), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_ifaces), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<FieldPointers>();

    m->n_nodes  = opp_get_uniform_local_size(g_m->n_nodes);
    m->n_cells  = opp_get_uniform_local_size(g_m->n_cells);
    m->n_ifaces = opp_get_uniform_local_size(g_m->n_ifaces);
    
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_ef      , m->c_ef      , g_m->n_cells , m->n_cells , DIM);
    opp_uniform_scatter_array(g_m->c_to_n    , m->c_to_n    , g_m->n_cells , m->n_cells , N_PER_C); 
    opp_uniform_scatter_array(g_m->c_to_c    , m->c_to_c    , g_m->n_cells , m->n_cells , NEIGHB_C); 
    opp_uniform_scatter_array(g_m->c_det     , m->c_det     , g_m->n_cells , m->n_cells , ALL_DET); 
    opp_uniform_scatter_array(g_m->c_vol     , m->c_vol     , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->c_sd      , m->c_sd      , g_m->n_cells , m->n_cells , N_PER_C*DIM); 
    opp_uniform_scatter_array(g_m->c_col     , m->c_col     , g_m->n_cells , m->n_cells , ONE);
    opp_uniform_scatter_array(g_m->c_id      , m->c_id      , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->n_bnd_pot , m->n_bnd_pot , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pot     , m->n_pot     , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_ion_den , m->n_ion_den , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pos     , m->n_pos     , g_m->n_nodes , m->n_nodes , DIM); 
    opp_uniform_scatter_array(g_m->n_vol     , m->n_vol     , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_type    , m->n_type    , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_color   , m->n_color   , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_bnd_pot , m->n_bnd_pot , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_id      , m->n_id      , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->if_to_c   , m->if_to_c   , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_to_n   , m->if_to_n   , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_v_norm , m->if_v_norm , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_u_norm , m->if_u_norm , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_norm   , m->if_norm   , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_area   , m->if_area   , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_dist   , m->if_dist   , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_n_pos  , m->if_n_pos  , g_m->n_ifaces, m->n_ifaces, DIM);       
    opp_uniform_scatter_array(g_m->if_id     , m->if_id     , g_m->n_ifaces, m->n_ifaces, ONE); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif
}

