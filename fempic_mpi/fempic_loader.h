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

#include "fempic_defs.h"
#include "fempic_ori/meshes.h"
#include "fempic_ori/particles.h"
// #include "fempic_ori/maths.h"

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
            if (c_centroid) delete[] c_centroid;
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
            c_centroid = nullptr;

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
            c_ef       = new double[n_cells * DIM];
            c_to_n     = new int[n_cells * N_PER_C];
            c_to_c     = new int[n_cells * NEIGHB_C];
            c_det      = new double[n_cells * ALL_DET]; // [alpha,beta,gamma,delta] * 4neighbours
            c_vol      = new double[n_cells];
            c_sd       = new double[n_cells * N_PER_C*DIM]; // arranged as [x,y,z] * 4 neighbours
            c_col      = new int[n_cells];
            c_id       = new int[n_cells];
            c_centroid = new double[n_cells * DIM];

            n_bnd_pot = new double[n_nodes];
            n_pot     = new double[n_nodes];
            n_ion_den = new double[n_nodes];
            n_pos     = new double[n_nodes * DIM];
            n_vol     = new double[n_nodes];
            n_type    = new int[n_nodes];
            n_color   = new int[n_nodes];
            n_id      = new int[n_nodes];

            if_to_c   = new int[n_ifaces];
            if_to_n   = new int[n_ifaces * N_PER_IF]; 
            if_v_norm = new double[n_ifaces * DIM]; 
            if_u_norm = new double[n_ifaces * DIM]; 
            if_norm   = new double[n_ifaces * DIM]; 
            if_area   = new double[n_ifaces]; 
            if_dist   = new int[n_ifaces];
            if_n_pos  = new double[n_ifaces * N_PER_IF * DIM];
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
        double *c_centroid = nullptr;

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
inline std::shared_ptr<FieldPointers> LoadMesh(std::string global_mesh, std::string inlet_mesh, 
                                        std::string wall_mesh, std::string cluster_type)
{
    std::shared_ptr<FieldPointers> mesh(new FieldPointers());

    std::shared_ptr<Volume> volume(new Volume());
    if (!LoadVolumeMesh(global_mesh, *(volume.get())) ||
        !(LoadSurfaceMesh(inlet_mesh, *(volume.get()),INLET, false)) ||
        !(LoadSurfaceMesh(wall_mesh, *(volume.get()), FIXED, false))) 
    {    
        return mesh;
    }

    if (OPP_rank == OPP_ROOT)
        volume->summarize(std::cout);

    mesh->n_nodes         = volume->nodes.size();
    mesh->n_cells         = volume->elements.size();
    mesh->n_ifaces        = volume->inlet_faces.size();

    printf("LoadMesh nodes %d cells %d ifaces %d\n", mesh->n_nodes, mesh->n_cells, mesh->n_ifaces);

    mesh->CreateMeshArrays();

    printf("LoadMesh CreateMeshArrays DONE\n");

    for (int n=0; n<mesh->n_nodes; n++)
    {
        mesh->n_pot[n] = 0.0f;
        mesh->n_ion_den[n] = 0.0f;

        Node &node = volume->nodes[n];
    
        for (int dim=0; dim<DIM; dim++)
            mesh->n_pos[n * DIM + dim] = node.pos[dim];
    
        mesh->n_vol[n] = node.volume;
        mesh->n_type[n]   = (int)node.type;
        mesh->n_id[n]     = n;

        if (n % 10000 == 0)
            printf("LoadMesh Mesh Nodes at %d\n", n);  
    }

    printf("LoadMesh Mesh Nodes Loaded\n");

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

        const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[0]].pos);
        const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[1]].pos);
        const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[2]].pos);
        const opp::Point3D *node3_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[3]].pos);
        const opp::Point3D c_pos = opp::getTetraCentroid3D(*node0_pos, *node1_pos, *node2_pos, *node3_pos);

        mesh->c_centroid[cID * DIM + 0] = c_pos.x;
        mesh->c_centroid[cID * DIM + 1] = c_pos.y;
        mesh->c_centroid[cID * DIM + 2] = c_pos.z;

        mesh->c_vol[cID] = tet.volume;
        mesh->c_id[cID] = cID;

        mesh->c_ef[cID * DIM + 0] = 0.0;
        mesh->c_ef[cID * DIM + 1] = 0.0;
        mesh->c_ef[cID * DIM + 2] = 0.0;        

        if (cID % 10000 == 0)
            printf("LoadMesh Mesh Cells at %d\n", cID);   
    }

    printf("LoadMesh Mesh Cells Loaded\n");

    for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        
        mesh->if_to_c[faceID] = face.cell_con;
        mesh->if_area[faceID] = face.area;
        mesh->if_dist[faceID] = 0;
        mesh->if_id[faceID] = faceID;

        for (int dim=0; dim<3; dim++)
        {
            mesh->if_v_norm[faceID * 3 + dim] = face.v[dim];   
            mesh->if_u_norm[faceID * 3 + dim] = face.u[dim];   
            mesh->if_norm[faceID * 3 + dim]   = face.normal[dim]; 
        }

        for (int n=0; n< N_PER_IF; n++)
        {
            mesh->if_to_n[faceID * 3 + n] = face.con[n]; 

            Node &node = volume->nodes[face.con[n]];
            for (int dim=0; dim<3; dim++)
            {
                mesh->if_n_pos[faceID * (N_PER_IF * DIM) + DIM * n + dim] = node.pos[dim];
            }
        }
    }

    printf("LoadMesh Mesh inlet faces Loaded\n");

#ifdef USE_MPI 

    // std::vector<std::pair<int, int>> face_pairs;    // face pair (to create a square using 2 triangle faces) vector
    // std::vector<opp::Point3D> face_centroids;            // centroids of the face pairs (should be the centroid of the square)
    // std::vector<int> cluster_assignments;
    // bool blockPartitioning = false;

    // if (cluster_type == "block") {
    //     opp_printf("INIT", "Using block cluster");
    //     blockPartitioning = true;

    //     std::map<int, std::vector<int>> node_face_con;

    //     // create a node to face connectivity mapping only for inlet face nodes
    //     for (int l=0; l < mesh->n_ifaces; l++) {
    //         Face &face = volume->inlet_faces[l];
    //         for (int v = 0; v < 3; v++)
    //             node_face_con[face.con[v]].emplace_back(l);
    //     }

    //     std::vector<bool> already_done(mesh->n_ifaces); // initializes with false by default

    //     for (int faceID = 0; faceID < mesh->n_ifaces; faceID++) {

    //         Face &face = volume->inlet_faces[faceID]; // the face of interest

    //         int v1,v2;
    //         for (int v = 0; v < 3; v++) { // iterate over the three vertices (nodes)
                
    //             if (already_done[faceID]) continue;

    //             switch(v) {
    //                 case 0: v1=1; v2=2; break;
    //                 case 1: v1=2; v2=0; break;
    //                 case 2: v1=0; v2=1; break;
    //             }

    //             for (int n = 0; n < 3; n++) { 

    //                 for (int m : node_face_con[face.con[n]]) {  // m are other face indices mapped with the node of the face of interest

    //                     if (faceID == m || already_done[m] || already_done[faceID]) continue;

    //                     Face &other = volume->inlet_faces[m];

    //                     bool matches[3] = { false, false, false };
    //                     int count = 0;
    //                     for (int k = 0; k < 3; k++) {
    //                         if (other.con[k] == face.con[v1] ||
    //                             other.con[k] == face.con[v2]) {
    //                                 count++;
    //                                 matches[k]=true;
    //                             }
    //                     }

    //                     if (count == 2) {

    //                         int nmatch_node = -1, nmatch_other_node = -1; // non-matching nodes of both the face triangles

    //                         nmatch_node = face.con[v];
    //                         for (int k = 0; k < 3; k++)
    //                             if(!matches[k]) 
    //                                 nmatch_other_node = other.con[k];
                            
    //                         // check whether the non-matching nodes are located diagonally (not aligning with x nor y)
    //                         if ((volume->nodes[nmatch_node].pos[0] != volume->nodes[nmatch_other_node].pos[0]) &&
    //                             (volume->nodes[nmatch_node].pos[1] != volume->nodes[nmatch_other_node].pos[1])) {
                                
    //                             face_pairs.push_back({faceID, m});
    //                             already_done[faceID] = true;
    //                             already_done[m] = true;

    //                             opp::Point3D centroid;
    //                             centroid.x = (volume->nodes[nmatch_node].pos[0] + volume->nodes[nmatch_other_node].pos[0]) / 2;
    //                             centroid.y = (volume->nodes[nmatch_node].pos[1] + volume->nodes[nmatch_other_node].pos[1]) / 2,
    //                             centroid.z = 0.0f;

    //                             face_centroids.push_back(centroid);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     node_face_con.clear();
    //     already_done.clear();

    //     cluster_assignments = opp::BlockCluster(face_centroids, OPP_comm_size);
    // }
    // else if (cluster_type == "k-means") {
    //     opp_printf("INIT", "Using k-means cluster");

    //     for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    //     {
    //         Face &face = volume->inlet_faces[faceID];
            
    //         const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[0]].pos);
    //         const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[1]].pos);
    //         const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[2]].pos);

    //         face_centroids.push_back(opp::getTriangleCentroid3D(*node0_pos, *node1_pos, *node2_pos));
    //     }

    //     cluster_assignments = opp::kMeansClustering3D(face_centroids, OPP_comm_size);
    // }
    // else {
    //     opp_abort("Error Cluster Type not defined");
    // }

    // bool color_found = false;

    // for (int cellID=0; cellID<mesh->n_cells; cellID++)
    // {
    //     Tetra &tet = volume->elements[cellID];
        
    //     const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[0]].pos);
    //     const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[1]].pos);
    //     const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[2]].pos);
    //     const opp::Point3D *node3_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[3]].pos);

    //     const opp::Point3D c_pos = opp::getTetraCentroid3D(*node0_pos, *node1_pos, *node2_pos, *node3_pos);

    //     color_found = false;

    //     for (int faceID=0; faceID<mesh->n_ifaces; faceID++) {
    //         Face &face = volume->inlet_faces[faceID];
            
    //         const opp::Point3D *face_node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[0]].pos);
    //         const opp::Point3D *face_node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[1]].pos);
    //         const opp::Point3D *face_node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[2]].pos);

    //         bool isInTriangle = opp::isPointInTriangle(c_pos, *face_node0_pos, *face_node1_pos, *face_node2_pos);

    //         if (isInTriangle) {
    //             int id = 0;
    //             if (blockPartitioning) {
    //                 for ( ; id < (int)face_pairs.size(); id++) {
    //                     if (face_pairs[id].first == faceID || face_pairs[id].second == faceID)
    //                         break;
    //                 }
    //             }
    //             else {
    //                 id = faceID;
    //             }

    //             mesh->c_col[cellID] = cluster_assignments[id];    
    //             color_found = true;
    //             break;
    //         }
    //     }

    //     if (!color_found) {
    //         opp_printf("Setup", "Error... Couldnt find colour for cell ,[%lf,%lf]", c_pos.x, c_pos.y);
    //     }   
    // }

    // // Print cluster assignments
    // if (false) // (OP_DEBUG)
    // {
    //     std::cout << "face_centroids : ";
    //     for (int i = 0; i < (int)face_centroids.size(); ++i) {
    //         std::cout << "[" << face_centroids[i].x << "," << face_centroids[i].y << "," << 
    //             face_centroids[i].z << "], ";
    //     }
    
    //     std::cout << std::endl << std::endl;

    //     std::cout << "cluster_assignments : ";
    //     for (int i = 0; i < (int)cluster_assignments.size(); ++i) {
    //         std::cout << cluster_assignments[i] << ",";
    //     }

    //     std::cout << std::endl << std::endl;
        
    //     std::vector<opp::Point3D> cluster_centroids = opp::calculateTriangleCentroids3D(face_centroids, 
    //                                                                         cluster_assignments);

    //     std::cout << "cluster_centroids : ";
    //     for (int i = 0; i < (int)cluster_centroids.size(); ++i) {
    //         std::cout << "[" << cluster_centroids[i].x << ","<< cluster_centroids[i].y << "," << 
    //             cluster_centroids[i].z << "], ";
    //     }

    //     std::cout << std::endl << std::endl;

    //     // DEBUG ONLY - to be removed start
    //         std::vector<double> face_areas_per_rank(OPP_comm_size, 0.0);
    //         for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    //         {
    //             Face &face = volume->inlet_faces[faceID];
                
    //             const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[0]].pos);
    //             const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[1]].pos);
    //             const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[face.con[2]].pos);

    //             face_areas_per_rank[cluster_assignments[faceID]] += 
    //                 std::abs(opp::calculateTriangleArea(*node0_pos, *node1_pos, *node2_pos) * 1000000);
    //         }
    //         double total = 0.0;
    //         std::cout << "face_areas_per_rank (um^2) -> ";
    //         for (int i = 0; i < OPP_comm_size; ++i) {
    //             total += face_areas_per_rank[i];
    //             std::cout << i << "|" << face_areas_per_rank[i] << " ";
    //         }
    //         std::cout << " = all|" << total << std::endl << std::endl;
    //     // DEBUG ONLY - to be removed end

    //     std::vector<int> c_counts(OPP_comm_size, 0);
    //     std::vector<int> face_counts(OPP_comm_size, 0);

    //     for (int n=0; n<mesh->n_cells; n++)
    //         c_counts[mesh->c_col[n]]++;

    //     for (int n=0; n<(int)cluster_assignments.size(); n++)
    //         face_counts[cluster_assignments[n]]++;

    //     std::string log = "";
    //     for (int n=0; n<OPP_comm_size; n++)
    //         log += std::to_string(n) + "|" + std::to_string(face_counts[n]) + "|" + 
    //             std::to_string(c_counts[n]) + " ";

    //     opp_printf("Setup", "Rank|faces|cells %s", log.c_str());
    // }   

#endif

    return mesh;
}

//*************************************************************************************************
inline void DistributeMeshOverRanks(std::shared_ptr<FieldPointers>& g_m, std::shared_ptr<FieldPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_nodes), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_ifaces), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<FieldPointers>();

    m->n_nodes  = opp_get_uniform_local_size(g_m->n_nodes);
    m->n_cells  = opp_get_uniform_local_size(g_m->n_cells);
    m->n_ifaces = opp_get_uniform_local_size(g_m->n_ifaces);
    
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_ef       , m->c_ef       , g_m->n_cells , m->n_cells , DIM);
    opp_uniform_scatter_array(g_m->c_to_n     , m->c_to_n     , g_m->n_cells , m->n_cells , N_PER_C); 
    opp_uniform_scatter_array(g_m->c_to_c     , m->c_to_c     , g_m->n_cells , m->n_cells , NEIGHB_C); 
    opp_uniform_scatter_array(g_m->c_det      , m->c_det      , g_m->n_cells , m->n_cells , ALL_DET); 
    opp_uniform_scatter_array(g_m->c_vol      , m->c_vol      , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->c_sd       , m->c_sd       , g_m->n_cells , m->n_cells , N_PER_C*DIM); 
    opp_uniform_scatter_array(g_m->c_col      , m->c_col      , g_m->n_cells , m->n_cells , ONE);
    opp_uniform_scatter_array(g_m->c_id       , m->c_id       , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->c_centroid , m->c_centroid , g_m->n_cells , m->n_cells , DIM); 
    opp_uniform_scatter_array(g_m->n_bnd_pot  , m->n_bnd_pot  , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pot      , m->n_pot      , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_ion_den  , m->n_ion_den  , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pos      , m->n_pos      , g_m->n_nodes , m->n_nodes , DIM); 
    opp_uniform_scatter_array(g_m->n_vol      , m->n_vol      , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_type     , m->n_type     , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_color    , m->n_color    , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_bnd_pot  , m->n_bnd_pot  , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_id       , m->n_id       , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->if_to_c    , m->if_to_c    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_to_n    , m->if_to_n    , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_v_norm  , m->if_v_norm  , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_u_norm  , m->if_u_norm  , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_norm    , m->if_norm    , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_area    , m->if_area    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_dist    , m->if_dist    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_n_pos   , m->if_n_pos   , g_m->n_ifaces, m->n_ifaces, N_PER_IF*DIM);       
    opp_uniform_scatter_array(g_m->if_id      , m->if_id      , g_m->n_ifaces, m->n_ifaces, ONE); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif
}

