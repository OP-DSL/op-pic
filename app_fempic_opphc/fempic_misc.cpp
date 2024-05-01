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
// USER WRITTEN CODE
//*********************************************

#define INJ_EXCESS 100
#define TOLERENCE 1E-12
#define PRECISION 12

//*************************************************************************************************
inline int init_inject_distributions(opp_dat if_dist_dat, opp_dat if_area_dat, 
                                        opp_dat dummy_random)
{
    if (OPP_DBG) opp_printf("init_inject_distributions", "START");

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

    if (OPP_DBG)
        opp_printf("init_inject_distributions", "RAND_FILE inj_count %d max_inj_count_per_face %d", 
            total_inject_count, max_inject_count_per_face);  
            
    // increase dummy random particle set size, to load the random numbers for particle injections
    opp_increase_particle_count(dummy_random->set, (total_inject_count + INJ_EXCESS));  

    int total_size = -1, fsize = -1, fdim = -1;
    FILE *fp = NULL;
    std::string rand_file_path = opp_params->get<OPP_STRING>("rand_file");

    // read from MPI ROOT and broadcast // alternatively could use HDF5 files
    if (OPP_rank == OPP_ROOT)
    {       
        if ((fp = fopen(rand_file_path.c_str(), "r")) == NULL)
        {
            opp_printf("init_inject_distributions", "Unable to open file %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2)
        {
            opp_printf("init_inject_distributions", "Error reading file data from %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        total_size = fsize * fdim;

        if (max_inject_count_per_face * dummy_random->dim > total_size)
        {
            opp_printf("init_inject_distributions", "dim and/or set_size issue in file %s\n", 
                rand_file_path.c_str());
            opp_abort();     
        }
    }

#ifdef USE_MPI
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
                opp_printf("init_inject_distributions", "Error reading from %s at index %d\n", 
                    rand_file_path.c_str(), n);
                opp_abort();
            }
        }

        fclose(fp);
    }

#ifdef USE_MPI
    // Load the whole file and bradcast
    // TODO : We can reduce communications by sending only the required size
    MPI_Bcast(dist, total_size, MPI_DOUBLE, OPP_ROOT, MPI_COMM_WORLD);
#endif

    if (if_area_dat->set->size > 0) 
    {
        double* random_dat = (double *)dummy_random->data;
        int* distribution  = (int *)if_dist_dat->data;
        int iface_index = 0, part_in_face_index = 0;

        // This trouble is only because we need mpi results to match with seq and others
        for (int i = 0; i < dummy_random->set->size; i++)
        {
            if (iface_index < if_dist_dat->set->size && i >= distribution[iface_index])
            {
                iface_index++; // check whether it is j or j-1
                part_in_face_index = 0; 
            }

            random_dat[i * 2 + 0] = dist[part_in_face_index * 2 + 0];
            random_dat[i * 2 + 1] = dist[part_in_face_index * 2 + 1];

            part_in_face_index++;
        }
    }
    
    delete[] dist;

    dummy_random->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    if_dist_dat->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    
    if (OPP_DBG) opp_printf("init_inject_distributions", "total_inject_count %d", total_inject_count);

    return total_inject_count;
}

//*************************************************************************************************
inline void print_per_cell_particle_counts(opp_dat c_part_count, opp_dat part_mesh_relation)     
{
    // validate the final particle counts in each cell
    opp_reset_dat(c_part_count, (char*)opp_zero_int16);
    for (int p = 0; p < part_mesh_relation->set->size; p++)
    {
        int c_index    = ((int *)part_mesh_relation->data)[p];
        ((int *)c_part_count->data)[c_index] += 1;
    }

#ifdef USE_MPI
    opp_mpi_print_dat_to_txtfile(c_part_count, "c_part_count.dat");
#else
    opp_print_dat_to_txtfile(c_part_count, "", "c_part_count.dat");
#endif
}

//*************************************************************************************************
inline std::string get_global_level_log(double max_n_charge_density, double max_n_potential, 
    int local_part_count, int local_parts_injected, int local_part_removed)
{
    std::string log = "";
    int64_t global_part_size = 0, global_inj_size = 0, global_removed = 0;
    int global_max_comm_iteration = 0;

#ifdef USE_MPI
    MPI_Reduce(&OPP_max_comm_iteration, &global_max_comm_iteration, 1, MPI_INT, MPI_MAX, OPP_ROOT, MPI_COMM_WORLD);

    int64_t temp_local_part_count     = (int64_t)local_part_count;
    int64_t temp_local_parts_injected = (int64_t)local_parts_injected;
    int64_t temp_local_part_removed   = (int64_t)local_part_removed;
    MPI_Reduce(&temp_local_part_count, &global_part_size, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&temp_local_parts_injected, &global_inj_size, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&temp_local_part_removed, &global_removed, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    
#else
    global_part_size = local_part_count;
    global_inj_size = local_parts_injected;
    global_removed = local_part_removed;
    global_max_comm_iteration = OPP_max_comm_iteration;
#endif

    log += std::string("\t np: ") + str(global_part_size, "%" PRId64);
    log += std::string(" (") + str(global_inj_size, "%" PRId64);
    log += std::string(" added, ") + str(global_removed, "%" PRId64);
    log += std::string(" removed)\t max den: ") + str(max_n_charge_density, "%2.25lE");
    log += std::string(" max |phi|: ") + str(max_n_potential, "%2.10lE");
    log += std::string(" max_comm_iteration: ") + str(global_max_comm_iteration, "%d");
    return log;
}

//*************************************************************************************************
inline int getGlobalSetSizes(opp_set set, std::vector<int>& counts_vec, std::vector<int>& ifaces_offsets) {

    counts_vec.clear();
    counts_vec.resize(OPP_comm_size, 0);

    ifaces_offsets.clear();
    ifaces_offsets.resize(OPP_comm_size, -1);

#ifdef USE_MPI
    MPI_Allgather(&(set->size), 1, MPI_INT, &counts_vec[0], 1, MPI_INT, MPI_COMM_WORLD);

    int count = 0;
    for (int i = 0; i < OPP_comm_size; i++) {
        ifaces_offsets[i] = count;
        count += counts_vec[i];
    }

    return count;
#else
    counts_vec[0] = set->size;
    ifaces_offsets[0] = 0;

    return set->size;
#endif
}

#ifdef USE_MPI
inline std::string convert_to_string(const double& double_value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(12) << double_value;
    return oss.str();
}

std::map<std::string, std::map<std::string, int>> get_strgrid_and_dims(const std::vector<opp::Point3D>& points, std::vector<int>& cell_counts) {

    std::map<std::string, std::map<std::string, int>> grid;

    for (int i = 0; i < (int)points.size(); i++) {
        grid[convert_to_string(points[i].y)][convert_to_string(points[i].x)] = i;
    }

    cell_counts[1] = (int)grid.size();
    for (const auto& outer_pair : grid) {
        const auto& inner_map = outer_pair.second;
        if ((int)inner_map.size() > cell_counts[0]) cell_counts[0] = (int)inner_map.size();
    }

    if ((int)points.size() != cell_counts[0] * cell_counts[1]) {
        opp_printf("Error", "Error expected_cell_count=%d x=%d y=%d", (int)points.size(), cell_counts[0], cell_counts[1]);
        opp_abort();
    }

    opp_printf("CENTROID_GRID", "SIZE x=%d y=%d", cell_counts[0], cell_counts[1]);
    // for (const auto& pair : grid) {
    //     std::string log = pair.first + " | ";
    //     for (const auto& pair_i : pair.second) log += pair_i.first + " , ";
    //     opp_printf("CENTROID_GRID", "%s", log.c_str());
    // }

    return grid;
}

inline void get_all_start_ends_mpi(const std::vector<int>& cell_counts, std::vector<int>& all_cell_starts, std::vector<int>& all_cell_ends) {

    const int ndim = 2;
    MPI_Comm comm_cart;
    int mpi_dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 1};
    int coords[3] = {0, 0, 0};
    int cell_starts[3] = {0, 0, 0}; // Holds the first cell this rank owns in each dimension.
    int cell_ends[3] = {1, 1, 1}; // Holds the last cell+1 this ranks owns in each dimension.

    MPI_Dims_create(OPP_comm_size, ndim, mpi_dims);

    std::vector<int> cell_count_ordering = reverse_argsort(cell_counts); // direction with most cells first to match mpi_dims order

    std::vector<int> mpi_dims_reordered(ndim);
    for (int dimx = 0; dimx < ndim; dimx++)  // reorder the mpi_dims to match the actual domain
        mpi_dims_reordered[cell_count_ordering[dimx]] = mpi_dims[dimx];

    MPI_Cart_create(MPI_COMM_WORLD, ndim, mpi_dims_reordered.data(), periods, 1, &comm_cart);
    MPI_Cart_get(comm_cart, ndim, mpi_dims, periods, coords);

    for (int dimx = 0; dimx < ndim; dimx++) 
        get_decomp_1d(mpi_dims[dimx], cell_counts[dimx], coords[dimx], &cell_starts[dimx], &cell_ends[dimx]);

    all_cell_starts.resize(OPP_comm_size * ndim);
    all_cell_ends.resize(OPP_comm_size * ndim);

    MPI_Allgather(cell_starts, ndim, MPI_INT, all_cell_starts.data(), ndim, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(cell_ends, ndim, MPI_INT, all_cell_ends.data(), ndim, MPI_INT, MPI_COMM_WORLD);

    // if (OPP_rank == OPP_ROOT)
    // {
    //     std::string log = "";
    //     for (int r = 0; r < OPP_comm_size; r++) {
    //         log += std::string("\nrank ") + std::to_string(r) + " start (";
    //         for (int d = 0; d < ndim; d++)
    //             log += std::to_string(all_cell_starts[r*ndim+d]) + ",";
    //         log += ") end (";
    //         for (int d = 0; d < ndim; d++)
    //             log += std::to_string(all_cell_ends[r*ndim+d]) + ",";
    //         log += ")";
    //     }
    //     opp_printf("get_all_start_ends_mpi", "%s", log.c_str());
    // }

    MPI_Comm_free(&comm_cart);
}

inline std::vector<int> fempicMPIDist(const std::map<std::string, std::map<std::string, int>>& grid, 
    std::vector<int>& cell_counts, const std::vector<int>& cell_starts, const std::vector<int>& cell_ends) {
    
    opp_printf("SETUP", "Distributing using fempicMPIDist");

    std::vector<int> assignments(cell_counts[0] * cell_counts[1]);

    int coord[2] = { -1, -1 };

    for (const auto& pair : grid) {
        coord[1]++;
        coord[0] = -1;
        
        for (const auto& pair_i : pair.second){
            
            coord[0]++;
            const int id = pair_i.second;
            bool assigned = false;

            for (int rank = 0; rank < OPP_comm_size; rank++) 
            {
                bool is_rank_suitable[2] = { true, true };

                for (int dimx = 0; dimx < 2; dimx++) 
                {
                    if ((cell_starts[2*rank+dimx] > coord[dimx]) || 
                        (cell_ends[2*rank+dimx] <= coord[dimx]))
                    {
                        is_rank_suitable[dimx] = false;
                        break;
                    }
                }

                if (is_rank_suitable[0] && is_rank_suitable[1])
                {
                    assignments[id] = rank;
                    assigned = true;
                    break;
                }
            }
            
            if (assigned == false) {
                opp_printf("Error", "Issue in assigning %s %s of id %d to grid", pair.first.c_str(), 
                                            pair_i.first.c_str(), pair_i.second);
            }
        }
    }
    
    return assignments;
}




//*************************************************************************************************
inline void genColoursForBlockPartition(opp_dat cell_colours, opp_dat cell_centroids, 
                        opp_dat iface_n_pos, opp_map iface_v_node_map) {
    
    opp_profiler->start("genColsForPart");

    std::vector<int> num_ifaces_vec, ifaces_offsets;
    const int g_iface_size = getGlobalSetSizes(iface_n_pos->set, num_ifaces_vec, ifaces_offsets);
    
    std::vector<int> g_face_rank_assignments(g_iface_size, -1);
    std::vector<double> g_if_npos_dat(g_iface_size * N_PER_IF * DIM);
    std::vector<MPI_Request> recv_requests;
    std::vector<int> g_if_to_n_map;

    std::vector<std::pair<int, int>> face_pairs;  // face pair (to create a square using 2 triangle faces) vector
    std::vector<opp::Point3D> face_centroids;     // centroids of the face pairs (should be the centroid of the square)
    std::map<std::string, std::map<std::string, int>> grid;
    std::vector<int> cell_counts = { 0, 0, 0 };

    // Gather iface_v_node_map and iface_n_pos to ROOT rank
    if (OPP_rank == OPP_ROOT) {
        
        int index = 0;
        recv_requests.resize(OPP_comm_size *2);
        g_if_to_n_map.resize(g_iface_size * N_PER_IF);

        for (int source_rank = 0; source_rank < OPP_comm_size; source_rank++) {
            
            int num_ifaces = num_ifaces_vec[source_rank];
     
            MPI_Irecv(&(g_if_to_n_map[index * N_PER_IF]), num_ifaces * N_PER_IF, MPI_INT, 
                        source_rank, 0, MPI_COMM_WORLD, &(recv_requests[source_rank*2]));
            MPI_Irecv(&(g_if_npos_dat[index * N_PER_IF * DIM]), num_ifaces * N_PER_IF * DIM, MPI_DOUBLE, 
                        source_rank, 1, MPI_COMM_WORLD, &(recv_requests[source_rank*2+1]));

            index += num_ifaces;
        }
    } 
    MPI_Send(iface_v_node_map->map, iface_v_node_map->from->size * N_PER_IF, MPI_INT, 0,
                0, MPI_COMM_WORLD);
    MPI_Send(iface_n_pos->data, iface_v_node_map->from->size * N_PER_IF * DIM, MPI_DOUBLE, 0,
                1, MPI_COMM_WORLD);
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUS_IGNORE);

    if (OPP_rank == OPP_ROOT) {

        // create a node to face connectivity mapping only for inlet face nodes
        std::map<int, std::vector<int>> node_face_con;
        for (int l=0; l < g_iface_size; l++) {
            for (int v = 0; v < N_PER_IF; v++) {
                const int n_index = g_if_to_n_map[l * N_PER_IF + v];
                node_face_con[n_index].emplace_back(l);
            }
        }

        std::vector<bool> already_done(g_iface_size); // initializes with false by default
        
        for (int faceID = 0; faceID < g_iface_size; faceID++) {

            int v1,v2;
            for (int v = 0; v < 3; v++) { // iterate over the three vertices (nodes)
                
                if (already_done[faceID]) continue;

                switch(v) {
                    case 0: v1=1; v2=2; break;
                    case 1: v1=2; v2=0; break;
                    case 2: v1=0; v2=1; break;
                }

                for (int n = 0; n < 3; n++) { 
                    
                    const int n_index = g_if_to_n_map[faceID * N_PER_IF + n];
                    for (int m : node_face_con[n_index]) {  // m are other face indices mapped with the node of the face of interest

                        if (faceID == m || already_done[m] || already_done[faceID]) continue;

                        bool matches[3] = { false, false, false };
                        int count = 0;
                        for (int k = 0; k < 3; k++) {
                            if (g_if_to_n_map[m * N_PER_IF + k] == g_if_to_n_map[faceID * N_PER_IF + v1] ||
                                g_if_to_n_map[m * N_PER_IF + k] == g_if_to_n_map[faceID * N_PER_IF + v2]) {
                                    count++;
                                    matches[k]=true;
                                }
                        }

                        if (count == 2) {

                            const int n_idx = faceID * N_PER_IF * DIM + v * DIM; // Matching node index in the if_npos dat
                            int othr_n_idx = -1;              // Non-matching node index in the if_npos dat

                            for (int k = 0; k < 3; k++)
                                if(!matches[k]) 
                                    othr_n_idx = m * N_PER_IF * DIM + k * DIM;
                            
                            // check whether the non-matching nodes are located diagonally (not aligning with x nor y)
                            if ((g_if_npos_dat[n_idx + 0] != g_if_npos_dat[othr_n_idx + 0]) &&
                                (g_if_npos_dat[n_idx + 1] != g_if_npos_dat[othr_n_idx + 1])) {
                                
                                face_pairs.push_back({faceID, m});
                                already_done[faceID] = true;
                                already_done[m] = true;

                                opp::Point3D centroid;
                                centroid.x = (g_if_npos_dat[n_idx + 0] + g_if_npos_dat[othr_n_idx + 0]) / 2;
                                centroid.y = (g_if_npos_dat[n_idx + 1] + g_if_npos_dat[othr_n_idx + 1]) / 2,
                                centroid.z = 0.0f;

                                face_centroids.push_back(centroid);
                            }
                        }
                    }
                }
            }
        }

        node_face_con.clear();
        already_done.clear();

        grid = get_strgrid_and_dims(face_centroids, cell_counts);
        
        for (int rank = 1; rank < OPP_comm_size; rank++) {
            MPI_Send(cell_counts.data(), cell_counts.size(), MPI_INT, rank, 1100, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(cell_counts.data(), cell_counts.size(), MPI_INT, 0, 1100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::vector<int> all_cell_starts, all_cell_ends;
    get_all_start_ends_mpi(cell_counts, all_cell_starts, all_cell_ends);

    if (OPP_rank == OPP_ROOT) {

        std::vector<int> cluster_assignments;
        std::string cluster_type = opp_params->get<OPP_STRING>("cluster");

        if (cluster_type == "block")
            cluster_assignments = opp::BlockCluster(face_centroids, OPP_comm_size);
        else if (cluster_type == "k-means")
            cluster_assignments = opp::kMeansClustering3D(face_centroids, OPP_comm_size);
        else if (cluster_type == "mpi-block")
            cluster_assignments = fempicMPIDist(grid, cell_counts, all_cell_starts, all_cell_ends);
        else
            opp_abort("Cluster type not recognized");

        for (size_t id = 0; id < cluster_assignments.size(); id++) {

            g_face_rank_assignments[face_pairs[id].first]  = cluster_assignments[id];
            g_face_rank_assignments[face_pairs[id].second] = cluster_assignments[id];
        }

        cluster_assignments.clear();
        face_pairs.clear();
    }

    // Distribute face node positions and face rank assignments over all ranks
    if (OPP_rank == OPP_ROOT) {

        for (int rank = 1; rank < OPP_comm_size; rank++) {

            MPI_Send(g_if_npos_dat.data(), g_iface_size * N_PER_IF * DIM, MPI_DOUBLE, rank,
                1000, MPI_COMM_WORLD);

            MPI_Send(g_face_rank_assignments.data(), g_iface_size, MPI_INT, rank,
                2000, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Status status1, status2;

        MPI_Recv(g_if_npos_dat.data(), g_iface_size * N_PER_IF * DIM, MPI_DOUBLE, 0, 1000, 
            MPI_COMM_WORLD, &status1);

        MPI_Recv(g_face_rank_assignments.data(), g_iface_size, MPI_INT, 0, 2000, 
            MPI_COMM_WORLD, &status2);
        
        int count;
        MPI_Get_count(&status1, MPI_CHAR, &count);  
        if ((size_t)count != g_iface_size * N_PER_IF * DIM * sizeof(double)) {
            opp_printf("ERROR", "g_if_npos_dat received count %d, expected %zu", 
                count, g_iface_size * N_PER_IF * DIM * sizeof(double));
            opp_abort();
        }
        MPI_Get_count(&status2, MPI_CHAR, &count);  
        if ((size_t)count != g_iface_size * sizeof(int)) {
            opp_printf("ERROR", "g_face_rank_assignments received count %d, expected %zu", 
                count, g_iface_size * sizeof(int));
            opp_abort();
        }    
    }

    // Assign the cell colours according to its centroid position relative to inlet faces
    bool color_found = false;
    int error_count = 0;
    for (int cellID=0; cellID<cell_centroids->set->size; cellID++)
    {
        const opp::Point3D *c_centroid = reinterpret_cast<opp::Point3D*>(&(((double*)cell_centroids->data)[cellID * DIM]));

        color_found = false;

        for (int faceID=0; faceID<g_iface_size; faceID++) {

            const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 0 * DIM]));
            const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 1 * DIM]));
            const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 2 * DIM]));

            if (opp::isPointInTriangle(*c_centroid, *node0_pos, *node1_pos, *node2_pos)) {

                ((int*)cell_colours->data)[cellID] = g_face_rank_assignments[faceID];    
                color_found = true;
                break;
            }
        }

        if (!color_found) {
            opp_printf("Setup", "Error... Couldnt find colour for cell ,[%2.25lE,%2.25lE]", c_centroid->x, c_centroid->y);
            error_count++;
        }   
    }

    opp_profiler->end("genColsForPart");

    if (OPP_DBG) opp_printf("Setup", "Cell Colouring Issue Count = %d", error_count);
}
#endif

// //*************************************************************************************************
// inline void genColoursForBlockPartition(opp_dat cell_colours, opp_dat cell_centroids, 
//                         opp_dat iface_n_pos, opp_map iface_v_node_map) {
    
//     opp_profiler->start("genColsForPart");

// #ifdef USE_MPI
//     std::vector<int> num_ifaces_vec, ifaces_offsets;
//     const int g_iface_size = getGlobalSetSizes(iface_n_pos->set, num_ifaces_vec, ifaces_offsets);
//     std::vector<int> g_face_rank_assignments(g_iface_size, -1);

//     std::vector<int> g_if_to_n_map;
//     std::vector<double> g_if_npos_dat(g_iface_size * N_PER_IF * DIM);
//     std::vector<MPI_Request> recv_requests;
//     std::map<int, std::vector<int>> node_face_con;

//     // Gather iface_v_node_map and iface_n_pos to ROOT rank
//     if (OPP_rank == OPP_ROOT) {
        
//         int index = 0;
//         recv_requests.resize(OPP_comm_size *2);
//         g_if_to_n_map.resize(g_iface_size * N_PER_IF);

//         for (int source_rank = 0; source_rank < OPP_comm_size; source_rank++) {
            
//             int num_ifaces = num_ifaces_vec[source_rank];
     
//             MPI_Irecv(&(g_if_to_n_map[index * N_PER_IF]), num_ifaces * N_PER_IF, MPI_INT, 
//                         source_rank, 0, MPI_COMM_WORLD, &(recv_requests[source_rank*2]));
//             MPI_Irecv(&(g_if_npos_dat[index * N_PER_IF * DIM]), num_ifaces * N_PER_IF * DIM, MPI_DOUBLE, 
//                         source_rank, 1, MPI_COMM_WORLD, &(recv_requests[source_rank*2+1]));

//             index += num_ifaces;
//         }
//     } 
//     MPI_Send(iface_v_node_map->map, iface_v_node_map->from->size * N_PER_IF, MPI_INT, 0,
//                 0, MPI_COMM_WORLD);
//     MPI_Send(iface_n_pos->data, iface_v_node_map->from->size * N_PER_IF * DIM, MPI_DOUBLE, 0,
//                 1, MPI_COMM_WORLD);
//     MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUS_IGNORE);

//     if (OPP_rank == OPP_ROOT) {

//         // create a node to face connectivity mapping only for inlet face nodes
//         node_face_con.clear();
//         for (int l=0; l < g_iface_size; l++) {
//             for (int v = 0; v < N_PER_IF; v++) {
//                 const int n_index = g_if_to_n_map[l * N_PER_IF + v];
//                 node_face_con[n_index].emplace_back(l);
//             }
//         }

//         std::vector<bool> already_done(g_iface_size); // initializes with false by default
//         std::vector<std::pair<int, int>> face_pairs;  // face pair (to create a square using 2 triangle faces) vector
//         std::vector<opp::Point3D> face_centroids;     // centroids of the face pairs (should be the centroid of the square)
        
//         for (int faceID = 0; faceID < g_iface_size; faceID++) {

//             int v1,v2;
//             for (int v = 0; v < 3; v++) { // iterate over the three vertices (nodes)
                
//                 if (already_done[faceID]) continue;

//                 switch(v) {
//                     case 0: v1=1; v2=2; break;
//                     case 1: v1=2; v2=0; break;
//                     case 2: v1=0; v2=1; break;
//                 }

//                 for (int n = 0; n < 3; n++) { 
                    
//                     const int n_index = g_if_to_n_map[faceID * N_PER_IF + n];
//                     for (int m : node_face_con[n_index]) {  // m are other face indices mapped with the node of the face of interest

//                         if (faceID == m || already_done[m] || already_done[faceID]) continue;

//                         bool matches[3] = { false, false, false };
//                         int count = 0;
//                         for (int k = 0; k < 3; k++) {
//                             if (g_if_to_n_map[m * N_PER_IF + k] == g_if_to_n_map[faceID * N_PER_IF + v1] ||
//                                 g_if_to_n_map[m * N_PER_IF + k] == g_if_to_n_map[faceID * N_PER_IF + v2]) {
//                                     count++;
//                                     matches[k]=true;
//                                 }
//                         }

//                         if (count == 2) {

//                             const int n_idx = faceID * N_PER_IF * DIM + v * DIM; // Matching node index in the if_npos dat
//                             int othr_n_idx = -1;              // Non-matching node index in the if_npos dat

//                             for (int k = 0; k < 3; k++)
//                                 if(!matches[k]) 
//                                     othr_n_idx = m * N_PER_IF * DIM + k * DIM;
                            
//                             // check whether the non-matching nodes are located diagonally (not aligning with x nor y)
//                             if ((g_if_npos_dat[n_idx + 0] != g_if_npos_dat[othr_n_idx + 0]) &&
//                                 (g_if_npos_dat[n_idx + 1] != g_if_npos_dat[othr_n_idx + 1])) {
                                
//                                 face_pairs.push_back({faceID, m});
//                                 already_done[faceID] = true;
//                                 already_done[m] = true;

//                                 opp::Point3D centroid;
//                                 centroid.x = (g_if_npos_dat[n_idx + 0] + g_if_npos_dat[othr_n_idx + 0]) / 2;
//                                 centroid.y = (g_if_npos_dat[n_idx + 1] + g_if_npos_dat[othr_n_idx + 1]) / 2,
//                                 centroid.z = 0.0f;

//                                 face_centroids.push_back(centroid);
//                             }
//                         }
//                     }
//                 }
//             }
//         }

//         node_face_con.clear();
//         already_done.clear();
        
//         std::vector<int> cluster_assignments = opp::BlockCluster(face_centroids, OPP_comm_size);
//         // std::vector<int> cluster_assignments = opp::kMeansClustering3D(face_centroids, OPP_comm_size);
//         // std::vector<int> cluster_assignments1 = fempicMPIDist(face_centroids, OPP_comm_size);

//         for (size_t id = 0; id < cluster_assignments.size(); id++) {

//             g_face_rank_assignments[face_pairs[id].first]  = cluster_assignments[id];
//             g_face_rank_assignments[face_pairs[id].second] = cluster_assignments[id];
//         }

//         cluster_assignments.clear();
//         face_pairs.clear();
//     }

//     // Distribute face node positions and face rank assignments over all ranks
//     if (OPP_rank == OPP_ROOT) {

//         for (int rank = 1; rank < OPP_comm_size; rank++) {

//             MPI_Send(g_if_npos_dat.data(), g_iface_size * N_PER_IF * DIM, MPI_DOUBLE, rank,
//                 1000, MPI_COMM_WORLD);

//             MPI_Send(g_face_rank_assignments.data(), g_iface_size, MPI_INT, rank,
//                 2000, MPI_COMM_WORLD);
//         }
//     }
//     else {
//         MPI_Status status1, status2;

//         MPI_Recv(g_if_npos_dat.data(), g_iface_size * N_PER_IF * DIM, MPI_DOUBLE, 0, 1000, 
//             MPI_COMM_WORLD, &status1);

//         MPI_Recv(g_face_rank_assignments.data(), g_iface_size, MPI_INT, 0, 2000, 
//             MPI_COMM_WORLD, &status2);
        
//         int count;
//         MPI_Get_count(&status1, MPI_CHAR, &count);  
//         if ((size_t)count != g_iface_size * N_PER_IF * DIM * sizeof(double)) {
//             opp_printf("ERROR", "g_if_npos_dat received count %d, expected %zu", 
//                 count, g_iface_size * N_PER_IF * DIM * sizeof(double));
//             opp_abort();
//         }
//         MPI_Get_count(&status2, MPI_CHAR, &count);  
//         if ((size_t)count != g_iface_size * sizeof(int)) {
//             opp_printf("ERROR", "g_face_rank_assignments received count %d, expected %zu", 
//                 count, g_iface_size * sizeof(int));
//             opp_abort();
//         }

//         // int error_code;
//         // MPI_Error_class(status1.MPI_ERROR, &error_code);
//         // if (error_code != MPI_SUCCESS) {
//         //     opp_printf("ERROR", "ERROR in g_if_npos_dat");
//         //     opp_abort();
//         // }
//         // MPI_Error_class(status2.MPI_ERROR, &error_code);
//         // if (error_code != MPI_SUCCESS) {
//         //     opp_printf("ERROR", "ERROR in g_face_rank_assignments");
//         //     opp_abort();
//         // }      
//     }

//     // Assign the cell colours according to its centroid position relative to inlet faces
//     bool color_found = false;
//     int error_count = 0;
//     for (int cellID=0; cellID<cell_centroids->set->size; cellID++)
//     {
//         const opp::Point3D *c_centroid = reinterpret_cast<opp::Point3D*>(&(((double*)cell_centroids->data)[cellID * DIM]));

//         color_found = false;

//         for (int faceID=0; faceID<g_iface_size; faceID++) {

//             const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 0 * DIM]));
//             const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 1 * DIM]));
//             const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(&(g_if_npos_dat[faceID * N_PER_IF * DIM + 2 * DIM]));

//             if (opp::isPointInTriangle(*c_centroid, *node0_pos, *node1_pos, *node2_pos)) {

//                 ((int*)cell_colours->data)[cellID] = g_face_rank_assignments[faceID];    
//                 color_found = true;
//                 break;
//             }
//         }

//         if (!color_found) {
//             opp_printf("Setup", "Error... Couldnt find colour for cell ,[%2.25lE,%2.25lE]", c_centroid->x, c_centroid->y);
//             error_count++;
//         }   
//     }

//     opp_profiler->end("genColsForPart");

//     if (OPP_DBG) opp_printf("Setup", "Cell Colouring Issue Count = %d", error_count);
// #endif
// }

//*************************************************************************************************
std::map<int, std::map<int,int>> particle_counts;

//*************************************************************************************************
inline void logSetSizeStatistics(opp_set set, int logBoundaryCount = 1) {

#ifdef USE_MPI    
    std::vector<int> count_per_iter(OPP_comm_size, 0);
    MPI_Gather(&(set->size), 1, MPI_INT, count_per_iter.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    int sum = 0, average = 0;

    auto& cc = particle_counts[OPP_main_loop_iter];
    for (int i = 0; i < OPP_comm_size; i++) {
        cc[count_per_iter[i]] = i;
        sum += count_per_iter[i];
    }
    average = sum / OPP_comm_size;

    std::string max_log = ""; 
    int counter_max = 1;
    for (auto it = cc.rbegin(); it != cc.rend(); ++it) {
        max_log += std::string(" ") + std::to_string(it->first) + "|" + std::to_string(it->second);
        if (counter_max == logBoundaryCount) break;
        counter_max++;
    }
    std::string min_log = "";
    int counter_min = 1;
    for (auto it = cc.begin(); it != cc.end(); ++it) {
        min_log += std::string(" ") + std::to_string(it->first) + "|" + std::to_string(it->second);
        if (counter_min == logBoundaryCount) break;
        counter_min++;
    }

    if (OPP_rank == 0) {
        opp_printf("ParticleSet", "sum %d average %d max [%s ] min [%s ]", 
            sum, average, max_log.c_str(), min_log.c_str());
    }
#endif
}

//*************************************************************************************************